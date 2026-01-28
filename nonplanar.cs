// Grasshopper Script Instance - NON-PLANAR EXTRUSION MULTIPLIER CALCULATION
// 
// INPUTS:
//   segments       - Tree  - Line segments in print order (any structure: branches, spirals, mixed)
//   nomGap         - Item  - Nominal gap between layers (default: 1.0mm)
//   extrusionWidth - Item  - Extrusion bead width in mm (default: 2.0mm)
//   samplesPerSeg  - Item  - Number of ray samples per segment (1-10). More samples smooth out
//                            kinks from non-planar polyline construction where Z varies along
//                            segments. Use 1 for planar layers, 3-5 for non-planar spirals.
//   minGapFraction - Item  - Minimum gap as fraction of nomGap to consider valid support (default: 0.4)
//                            Filters out same-layer hits in spirals. Use 0 to disable.
//   minMult        - Item  - Minimum multiplier clamp (default: 0.003)
//   maxMult        - Item  - Maximum multiplier clamp (default: 3.0)
//   showDebug      - Item  - Show debug visualization (default: false)
//
// OUTPUTS:
//   multipliers   - Multipliers for each segment (matches input tree structure)
//   summary       - Statistics report
//   debugLines    - Debug: vertical lines showing gap from segment midpoint to support
//   beadMesh      - Debug: the mesh used for ray intersection
//
// ALGORITHM:
// Mesh-based ray casting approach:
//
// Pass 1 - Mesh Construction:
//   Build a single mesh where each segment becomes a quad face (ribbon).
//   Face index equals segment index for trivial temporal ordering.
//   Quads are extended slightly to ensure overlap at curves.
//
// Pass 2 - Ray Casting:
//   For each segment, cast rays at samplesPerSeg positions along the segment.
//   At each position, cast 3 rays: center and ±1/3 extrusionWidth perpendicular to segment.
//   This captures gap variation across bead width on angled walls.
//   Filter hits to require minimum Z gap (filters same-layer hits in spirals).
//   Use actual mesh hit Z for accurate gap calculation with extended quads.
//   Average all valid gaps across samples and offsets.
//   Multiplier = avgGap / nomGap
//
// VERSION HISTORY:
//   v1 - Original mesh ray cast approach
//   v2 - Fixed tilted quad Z interpolation (was using mesh hit Z instead of segment Z)
//   v3 - Added minGapFraction to filter same-layer hits in spirals; consolidated debug output
//   v4 - Added perpendicular offset rays at ±1/3 width to handle angled walls and narrow mesh gaps
//   v5 - Fixed extended quad Z bug: use actual mesh hit Z instead of interpolating segment Z
//        (interpolation was wrong for hits on extended portions of mesh quads)
//   v6 - Fixed mesh winding order for zig-zag paths: ensure consistent perpendicular direction
//        so all faces have upward normals (rays were missing faces with downward normals)

#region Usings
using System;
using System.Linq;
using System.Collections.Generic;
using System.Diagnostics;
using System.Text;
using System.Threading;
using System.Threading.Tasks;

using Rhino.Geometry;
using Rhino.Geometry.Intersect;

using Grasshopper;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Data;
using Grasshopper.Kernel.Types;
#endregion

public class Script_Instance : GH_ScriptInstance
{
    private void RunScript(
		DataTree<Line> segments,
		double nomGap,
		double extrusionWidth,
		int samplesPerSeg,
		double minGapFraction,
		double minMult,
		double maxMult,
		bool showDebug,
		ref object multipliers,
		ref object summary,
		ref object debugLines,
		ref object beadMesh)
    {
        // Apply defaults
        if (nomGap <= 0) nomGap = 1.0;
        if (extrusionWidth <= 0) extrusionWidth = 2.0;
        if (samplesPerSeg < 1) samplesPerSeg = 1;
        if (samplesPerSeg > 10) samplesPerSeg = 10;
        if (minGapFraction < 0) minGapFraction = 0;
        if (minGapFraction > 0.9) minGapFraction = 0.9;
        // Default minGapFraction to 0.4 if not specified (input is 0)
        if (minGapFraction == 0) minGapFraction = 0.4;
        if (minMult <= 0) minMult = 0.003;
        if (maxMult <= 0) maxMult = 3.0;
        
        // Minimum Z gap to consider valid support - filters same-layer hits in spirals
        double minZGap = nomGap * minGapFraction;
        
        // Perpendicular offsets: center and ±1/3 of extrusion width
        double perpOffset = extrusionWidth / 3.0;
        double[] perpOffsets = new double[] { -perpOffset, 0, perpOffset };
        
        var timer = Stopwatch.StartNew();
        
        // Flatten all segments into global list with branch info for output reconstruction
        IList<GH_Path> paths = segments.Paths;
        int pathCount = paths.Count;
        
        var allSegments = new List<Line>();
        var branchStarts = new int[pathCount];
        var branchCounts = new int[pathCount];
        
        for (int b = 0; b < pathCount; b++)
        {
            branchStarts[b] = allSegments.Count;
            List<Line> branch = segments.Branch(paths[b]);
            branchCounts[b] = branch.Count;
            allSegments.AddRange(branch);
        }
        
        int totalSegments = allSegments.Count;
        
        // Handle empty input
        if (totalSegments == 0)
        {
            multipliers = new DataTree<double>();
            summary = "No segments provided.";
            debugLines = new DataTree<Line>();
            beadMesh = null;
            return;
        }
        
        // ============================================================
        // PASS 1: Build mesh where each face represents a segment bead
        //         Also compute perpendicular unit vectors for ray offsets
        // ============================================================
        double halfWidth = extrusionWidth / 2.0;
        
        // Pre-allocate mesh arrays
        var vertices = new Point3d[totalSegments * 4];
        var faces = new MeshFace[totalSegments];
        
        // Store perpendicular unit vectors for each segment (needed in Pass 2)
        var perpUnits = new Vector3d[totalSegments];
        
        // Extend each quad slightly beyond segment endpoints to ensure overlap at curves
        double extendLength = extrusionWidth * 0.15;
        
        Parallel.For(0, totalSegments, segIdx =>
        {
            Line seg = allSegments[segIdx];
            Point3d p0 = seg.From;
            Point3d p1 = seg.To;
            
            // Direction along segment (full 3D for extension)
            Vector3d along = new Vector3d(p1.X - p0.X, p1.Y - p0.Y, p1.Z - p0.Z);
            double len = along.Length;
            
            Vector3d alongUnit;
            Vector3d perpUnit;
            if (len > 1e-10)
            {
                alongUnit = along / len;
                // Perpendicular in XY plane (for width)
                Vector3d alongXY = new Vector3d(along.X, along.Y, 0);
                double lenXY = alongXY.Length;
                if (lenXY > 1e-10)
                {
                    alongXY = alongXY / lenXY;
                    perpUnit = new Vector3d(-alongXY.Y, alongXY.X, 0);
                }
                else
                {
                    // Vertical segment
                    perpUnit = new Vector3d(1, 0, 0);
                }
            }
            else
            {
                // Degenerate segment - make a small square
                alongUnit = new Vector3d(1, 0, 0);
                perpUnit = new Vector3d(0, 1, 0);
            }
            
            // Ensure consistent perpendicular direction for mesh winding order
            // All faces should have normals pointing UP for reliable ray intersection
            // Flip perpUnit if it points in negative Y (or negative X when Y is zero)
            if (perpUnit.Y < 0 || (perpUnit.Y == 0 && perpUnit.X < 0))
            {
                perpUnit = -perpUnit;
            }

            // Store perpendicular unit vector for ray casting
            perpUnits[segIdx] = perpUnit;
            
            // Perpendicular scaled to half width for mesh quad
            Vector3d perp = perpUnit * halfWidth;
            
            // Extend endpoints along segment direction
            Point3d p0ext = p0 - alongUnit * extendLength;
            Point3d p1ext = p1 + alongUnit * extendLength;
            
            // Create quad vertices following segment's actual Z at each end
            // Winding order: counter-clockwise when viewed from above (normal points up)
            int vi = segIdx * 4;
            vertices[vi + 0] = new Point3d(p0ext.X - perp.X, p0ext.Y - perp.Y, p0ext.Z);
            vertices[vi + 1] = new Point3d(p0ext.X + perp.X, p0ext.Y + perp.Y, p0ext.Z);
            vertices[vi + 2] = new Point3d(p1ext.X + perp.X, p1ext.Y + perp.Y, p1ext.Z);
            vertices[vi + 3] = new Point3d(p1ext.X - perp.X, p1ext.Y - perp.Y, p1ext.Z);
            
            faces[segIdx] = new MeshFace(vi, vi + 1, vi + 2, vi + 3);
        });
        
        // Assemble mesh
        var mesh = new Mesh();
        mesh.Vertices.AddVertices(vertices);
        mesh.Faces.AddFaces(faces);
        mesh.Normals.ComputeNormals();
        mesh.Compact();
        
        var meshBuildTime = timer.ElapsedMilliseconds;
        
        // ============================================================
        // PASS 2: Ray cast from each segment to find support below
        // ============================================================
        var allMults = new double[totalSegments];
        var allDebugLines = showDebug ? new Line[totalSegments] : null;
        
        // Statistics
        int withSupport = 0, withoutSupport = 0, skippedSameLayer = 0;
        double minGapFound = double.MaxValue, maxGapFound = double.MinValue;
        double minMultFound = double.MaxValue, maxMultFound = double.MinValue;
        double multSum = 0;
        int clampedCount = 0;
        object statsLock = new object();
        
        // Find Z range for ray length
        double minZ = double.MaxValue;
        foreach (var seg in allSegments)
        {
            if (seg.From.Z < minZ) minZ = seg.From.Z;
            if (seg.To.Z < minZ) minZ = seg.To.Z;
        }
        
        // Pre-compute sample positions along segment (evenly spaced, avoiding endpoints)
        double[] sampleTs = new double[samplesPerSeg];
        for (int i = 0; i < samplesPerSeg; i++)
        {
            // Space samples evenly from 0.1 to 0.9 (avoiding exact endpoints)
            sampleTs[i] = 0.1 + 0.8 * i / Math.Max(1, samplesPerSeg - 1);
        }
        if (samplesPerSeg == 1) sampleTs[0] = 0.5; // Single sample at midpoint
        
        int totalRaysPerSegment = samplesPerSeg * perpOffsets.Length;
        
        Parallel.For(0, totalSegments, segIdx =>
        {
            Line seg = allSegments[segIdx];
            Vector3d perpUnit = perpUnits[segIdx];
            
            // Sample at multiple points along segment and across bead width
            int validSamples = 0;
            double gapSum = 0;
            double bestGap = double.MaxValue;
            int localSkippedSameLayer = 0;
            
            foreach (double t in sampleTs)
            {
                Point3d centerPt = seg.PointAt(t);
                
                // Cast rays at center and perpendicular offsets
                foreach (double offset in perpOffsets)
                {
                    Point3d samplePt = new Point3d(
                        centerPt.X + perpUnit.X * offset,
                        centerPt.Y + perpUnit.Y * offset,
                        centerPt.Z
                    );
                    
                    // Cast ray downward from sample point
                    Point3d rayEnd = new Point3d(samplePt.X, samplePt.Y, minZ - 10);
                    Line ray = new Line(samplePt, rayEnd);
                    
                    // Find all mesh intersections
                    Point3d[] hits;
                    int[] faceIds;
                    hits = Intersection.MeshLine(mesh, ray, out faceIds);
                    
                    // Find highest valid hit for this sample
                    double sampleBestZ = double.MinValue;
                    
                    if (hits != null)
                    {
                        // Find highest hit from any earlier segment that's far enough below
                        for (int i = 0; i < hits.Length; i++)
                        {
                            int hitFaceIdx = faceIds[i];
                            
                            // Must be from earlier segment
                            if (hitFaceIdx < segIdx)
                            {
                                // Use actual mesh hit Z (not interpolated segment Z)
                                // This correctly handles extended mesh quads where Z extends
                                // beyond the original segment endpoints
                                double hitZ = hits[i].Z;

                                // Must be far enough below us (filters same-layer hits in spirals)
                                if (hitZ < samplePt.Z - minZGap)
                                {
                                    if (hitZ > sampleBestZ)
                                    {
                                        sampleBestZ = hitZ;
                                    }
                                }
                                else if (hitZ < samplePt.Z)
                                {
                                    // This hit was below us but filtered out by minZGap
                                    localSkippedSameLayer++;
                                }
                            }
                        }
                    }
                    
                    if (sampleBestZ > double.MinValue)
                    {
                        double sampleGap = samplePt.Z - sampleBestZ;
                        gapSum += sampleGap;
                        validSamples++;
                        
                        // Track smallest gap for debug visualization
                        if (sampleGap < bestGap)
                        {
                            bestGap = sampleGap;
                        }
                    }
                }
            }
            
            // Use midpoint for debug line
            Point3d midPt = seg.PointAt(0.5);
            
            double mult;
            if (validSamples > 0)
            {
                double avgGap = gapSum / validSamples;
                double rawMult = avgGap / nomGap;
                mult = Math.Max(minMult, Math.Min(maxMult, rawMult));
                
                if (Math.Abs(mult - rawMult) > 0.1)
                    Interlocked.Increment(ref clampedCount);
                
                Interlocked.Increment(ref withSupport);
                lock(statsLock)
                {
                    if (avgGap < minGapFound) minGapFound = avgGap;
                    if (avgGap > maxGapFound) maxGapFound = avgGap;
                    if (mult < minMultFound) minMultFound = mult;
                    if (mult > maxMultFound) maxMultFound = mult;
                    multSum += mult;
                    skippedSameLayer += localSkippedSameLayer;
                }
                
                if (showDebug)
                {
                    // Show line from midpoint down by the calculated gap
                    double displayGap = mult * nomGap;
                    allDebugLines[segIdx] = new Line(midPt, new Point3d(midPt.X, midPt.Y, midPt.Z - displayGap));
                }
            }
            else
            {
                mult = 1.0;
                Interlocked.Increment(ref withoutSupport);
                lock(statsLock)
                {
                    if (mult < minMultFound) minMultFound = mult;
                    if (mult > maxMultFound) maxMultFound = mult;
                    multSum += mult;
                    skippedSameLayer += localSkippedSameLayer;
                }
                
                if (showDebug)
                {
                    // No support found - show nominal gap
                    allDebugLines[segIdx] = new Line(midPt, new Point3d(midPt.X, midPt.Y, midPt.Z - nomGap));
                }
            }
            
            allMults[segIdx] = mult;
        });
        
        timer.Stop();
        
        // ============================================================
        // Reconstruct output trees to match input branch structure
        // ============================================================
        var multiplierTree = new DataTree<double>();
        var debugLineTree = new DataTree<Line>();
        
        for (int b = 0; b < pathCount; b++)
        {
            GH_Path path = paths[b];
            int start = branchStarts[b];
            int count = branchCounts[b];
            
            // Extract multipliers for this branch
            var branchMults = new double[count];
            Array.Copy(allMults, start, branchMults, 0, count);
            multiplierTree.AddRange(branchMults, path);
            
            if (showDebug)
            {
                var lines = new Line[count];
                Array.Copy(allDebugLines, start, lines, 0, count);
                debugLineTree.AddRange(lines, path);
            }
        }
        
        // Generate summary
        var report = new StringBuilder();
        report.AppendLine("=== NON-PLANAR EXTRUSION MULTIPLIER (MESH RAY CAST v6) ===");
        report.AppendLine($"Segments: {totalSegments:N0} | Mesh faces: {mesh.Faces.Count:N0}");
        report.AppendLine($"Rays per segment: {totalRaysPerSegment} ({samplesPerSeg} along × 3 across) | Total rays: {totalSegments * totalRaysPerSegment:N0}");
        report.AppendLine($"Time: {timer.ElapsedMilliseconds}ms (mesh build: {meshBuildTime}ms)");
        report.AppendLine($"Settings: nomGap={nomGap}mm, extrusionWidth={extrusionWidth}mm");
        report.AppendLine($"Perpendicular offsets: ±{perpOffset:F2}mm (1/3 of extrusion width)");
        report.AppendLine($"Min gap filter: {minGapFraction:P0} of nomGap = {minZGap:F3}mm");
        report.AppendLine($"Multiplier clamp: [{minMult:F3}, {maxMult}]");
        report.AppendLine();
        
        report.AppendLine($"With support: {withSupport:N0} ({100.0*withSupport/totalSegments:F1}%)");
        report.AppendLine($"Without support (first layer): {withoutSupport:N0}");
        if (skippedSameLayer > 0)
            report.AppendLine($"Same-layer hits filtered: {skippedSameLayer:N0}");
        
        if (withSupport > 0)
        {
            report.AppendLine($"Gap range: {minGapFound:F3} - {maxGapFound:F2}mm");
            report.AppendLine($"Multiplier range: {minMultFound:F3} - {maxMultFound:F2}");
            report.AppendLine($"Multiplier avg: {multSum/totalSegments:F3}");
            if (clampedCount > 0)
                report.AppendLine($"Clamped: {clampedCount:N0}");
        }
        
        // Set outputs
        multipliers = multiplierTree;
        summary = report.ToString();
        debugLines = debugLineTree;
        beadMesh = showDebug ? mesh : null;
    }
}