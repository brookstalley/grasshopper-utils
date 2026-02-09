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
//   baseZ          - Item  - Z height of base-to-form transition. Segments above this Z whose
//                            closest support is below this Z use default multiplier (1.0) for
//                            strong adhesion. 0 = disabled. (default: 0)
//
// OUTPUTS:
//   multipliers   - Multipliers for each segment (matches input tree structure)
//   summary       - Statistics report
//   debugLines    - Debug: vertical lines showing gap from segment midpoint to support
//   beadMesh      - Debug: the mesh used for ray intersection
//
// GRASSHOPPER WIRING (component nickname: "Nonplanar"):
//   Inputs:
//     segments       ← "Print segments" Curve parameter (tree of line-like curves)
//     nomGap         ← Panel (value: "1")
//     extrusionWidth ← Panel (value: "2")
//     samplesPerSeg  ← Number Slider (integer, value: 9)
//     minGapFraction ← Number Slider (value: 0.1)
//     minMult        ← "Min extrusion mult" Number Slider (value: 0.25)
//     maxMult        ← "Max extrusion mult" Number Slider (value: 1.75)
//     showDebug      ← Boolean Toggle (value: True)
//     baseZ          ← Panel (value: "0", or Z height of base top)
//   Outputs:
//     multipliers    → Number parameter, Bounds component, Gradient component, Panel
//     summary        → "nonplanar-summary" Panel
//     debugLines     → Line parameter
//     beadMesh       → Mesh parameter
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
//   At each position, cast center ray first. If center misses, cast ±20% extrusionWidth
//   perpendicular offset rays as fallback and use closest support found.
//   Filter hits to require minimum Z gap (filters same-layer hits in spirals).
//   Use actual mesh hit Z for accurate gap calculation with extended quads.
//   Average valid gaps across sample positions along segment.
//   Multiplier = avgGap / nomGap
//   Overhangs (rawMult > maxMult) and cross-section transitions (baseAdhesion) use default 1.0
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
//   v7 - Fixed edge blowout: perpendicular offset rays at form edges would miss nearby support
//        and find distant faces far below, dragging the average wildly wrong. Now center ray
//        is primary; offset rays (reduced to ±20% width) only cast as fallback when center
//        misses. Take closest support from fallback rays instead of averaging all rays.
//        Overhangs (rawMult > maxMult) use default multiplier 1.0 instead of clamping to max.
//        Added baseZ: Z-height threshold for base adhesion. Segments above baseZ supported
//        by faces below baseZ use default multiplier for strong adhesion at transitions.

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
		DataTree<object> segments,
		double nomGap,
		double extrusionWidth,
		int samplesPerSeg,
		double minGapFraction,
		double minMult,
		double maxMult,
		bool showDebug,
		object baseZInput,
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

        // Convert baseZ from object (Generic Data parameter)
        double baseZ = 0;
        if (baseZInput is double) baseZ = (double)baseZInput;
        else if (baseZInput is int) baseZ = (int)baseZInput;
        else if (baseZInput is GH_Number) baseZ = ((GH_Number)baseZInput).Value;
        else if (baseZInput != null) double.TryParse(baseZInput.ToString(), out baseZ);

        // Minimum Z gap to consider valid support - filters same-layer hits in spirals
        double minZGap = nomGap * minGapFraction;
        
        // Perpendicular offset for fallback rays: ±20% of extrusion width
        double perpOffset = extrusionWidth * 0.2;
        
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
            var branch = segments.Branch(paths[b]);
            for (int i = 0; i < branch.Count; i++)
            {
                object item = branch[i];
                if (item is Line)
                    allSegments.Add((Line)item);
                else if (item is GH_Line)
                    allSegments.Add(((GH_Line)item).Value);
                else if (item is Rhino.Geometry.LineCurve)
                    allSegments.Add(((Rhino.Geometry.LineCurve)item).Line);
            }
            branchCounts[b] = allSegments.Count - branchStarts[b];
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
        int withSupport = 0, withoutSupport = 0, skippedSameLayer = 0, invalidSamples = 0, adhesionOverride = 0;
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

        Parallel.For(0, totalSegments, segIdx =>
        {
            Line seg = allSegments[segIdx];
            Vector3d perpUnit = perpUnits[segIdx];
            
            // Sample at multiple points along segment; center ray first with offset fallback
            int validSamples = 0;
            double gapSum = 0;
            double bestGap = double.MaxValue;
            int localSkippedSameLayer = 0;
            int localInvalidSamples = 0;
            int localCrossZ = 0;

            // Rays: center first, then offset rays as fallback if center misses
            double[] rayOffsets = new double[] { 0, -perpOffset, perpOffset };

            foreach (double t in sampleTs)
            {
                Point3d centerPt = seg.PointAt(t);
                double posGap = double.MaxValue;
                bool posHasHit = false;

                for (int ri = 0; ri < rayOffsets.Length; ri++)
                {
                    // Skip offset rays if center ray (ri==0) already found support
                    if (ri == 1 && posHasHit) break;

                    Point3d samplePt = new Point3d(
                        centerPt.X + perpUnit.X * rayOffsets[ri],
                        centerPt.Y + perpUnit.Y * rayOffsets[ri],
                        centerPt.Z);

                    // Cast ray downward from sample point
                    Point3d rayEnd = new Point3d(samplePt.X, samplePt.Y, minZ - 10);
                    Line ray = new Line(samplePt, rayEnd);

                    // Find all mesh intersections
                    Point3d[] hits;
                    int[] faceIds;
                    hits = Intersection.MeshLine(mesh, ray, out faceIds);

                    if (hits != null)
                    {
                        for (int i = 0; i < hits.Length; i++)
                        {
                            int hitFaceIdx = faceIds[i];

                            // Must be from earlier segment
                            if (hitFaceIdx < segIdx)
                            {
                                double hitZ = hits[i].Z;

                                // Must be far enough below (filters same-layer hits)
                                if (hitZ < samplePt.Z - minZGap)
                                {
                                    double g = samplePt.Z - hitZ;
                                    if (g < posGap)
                                    {
                                        posGap = g;
                                        posHasHit = true;
                                    }
                                }
                                else if (hitZ < samplePt.Z)
                                {
                                    localSkippedSameLayer++;
                                }
                            }
                        }
                    }
                }

                if (posHasHit)
                {
                    gapSum += posGap;
                    validSamples++;
                    if (posGap < bestGap) bestGap = posGap;
                    // Cross-Z-boundary: segment above baseZ, support below baseZ
                    if (baseZ > 0 && centerPt.Z > baseZ && (centerPt.Z - posGap) < baseZ)
                        localCrossZ++;
                }
                else
                {
                    localInvalidSamples++;
                }
            }
            
            // Use midpoint for debug line
            Point3d midPt = seg.PointAt(0.5);
            
            double mult;
            if (validSamples > 0)
            {
                double avgGap = gapSum / validSamples;
                double rawMult = avgGap / nomGap;
                if (baseZ > 0 && localCrossZ * 2 > validSamples)
                {
                    mult = 1.0; // Base adhesion: support is below baseZ threshold
                    Interlocked.Increment(ref adhesionOverride);
                }
                else if (rawMult > maxMult)
                    mult = 1.0; // Overhang: gap too large, use default flow
                else
                    mult = Math.Max(minMult, rawMult);

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
                    invalidSamples += localInvalidSamples;
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
                    // Don't count invalid samples for no-support segments;
                    // they're already reported as "without support"
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
        report.AppendLine("=== NON-PLANAR EXTRUSION MULTIPLIER (MESH RAY CAST v7) ===");

        report.AppendLine(String.Format("Segments: {0:N0} | Mesh faces: {1:N0}", totalSegments, mesh.Faces.Count));
        report.AppendLine(String.Format("Samples: {0} per segment (center ray + +/-{1:F2}mm fallback offsets)", samplesPerSeg, perpOffset));
        report.AppendLine(String.Format("Time: {0}ms (mesh build: {1}ms)", timer.ElapsedMilliseconds, meshBuildTime));
        report.AppendLine(String.Format("Settings: nomGap={0}mm, extrusionWidth={1}mm", nomGap, extrusionWidth));
        report.AppendLine(String.Format("Min gap filter: {0:P0} of nomGap = {1:F3}mm", minGapFraction, minZGap));
        report.AppendLine(String.Format("Multiplier clamp: [{0:F3}, {1}] (overhang > max = 1.0)", minMult, maxMult));
        if (baseZ > 0)
            report.AppendLine(String.Format("Base adhesion: Z < {0}mm = 1.0", baseZ));
        report.AppendLine();

        report.AppendLine(String.Format("With support: {0:N0} ({1:F1}%)", withSupport, 100.0*withSupport/totalSegments));
        report.AppendLine(String.Format("Without support (first layer): {0:N0}", withoutSupport));
        if (invalidSamples > 0)
            report.AppendLine(String.Format("Invalid samples (no ray hit): {0:N0}", invalidSamples));
        if (skippedSameLayer > 0)
            report.AppendLine(String.Format("Same-layer hits filtered: {0:N0}", skippedSameLayer));

        if (withSupport > 0)
        {
            report.AppendLine(String.Format("Gap range: {0:F3} - {1:F2}mm", minGapFound, maxGapFound));
            report.AppendLine(String.Format("Multiplier range: {0:F3} - {1:F2}", minMultFound, maxMultFound));
            report.AppendLine(String.Format("Multiplier avg: {0:F3}", multSum/totalSegments));
            if (clampedCount > 0)
                report.AppendLine(String.Format("Clamped: {0:N0}", clampedCount));
            if (adhesionOverride > 0)
                report.AppendLine(String.Format("Base adhesion override: {0:N0}", adhesionOverride));
        }
        
        // Set outputs
        multipliers = multiplierTree;
        summary = report.ToString();
        debugLines = debugLineTree;
        beadMesh = showDebug ? mesh : null;
    }
}