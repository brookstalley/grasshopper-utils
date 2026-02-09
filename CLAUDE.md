# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This repository contains C# script components for Grasshopper (Rhino 3D's visual programming environment). These are not standalone applications—they run inside Grasshopper's embedded C# script editor.

## Build & Development

**No external build system.** Scripts are edited and executed directly in Grasshopper's C# Script component. To test changes:
1. Open Rhino 3D with Grasshopper
2. Create or open a C# Script component
3. Paste the script content
4. Connect inputs and run

## Architecture

### nonplanar.cs - Extrusion Multiplier Calculator

Calculates extrusion flow multipliers for non-planar 3D printing by measuring gaps between layers using ray casting.

**Two-Pass Algorithm:**

1. **Pass 1 - Mesh Construction** (parallel): Each segment becomes a quad face (ribbon mesh). Face index equals segment index for temporal ordering. Quads extend slightly beyond endpoints for overlap at curves.

2. **Pass 2 - Ray Casting** (parallel): For each segment, cast rays downward at multiple sample points. At each position, cast 3 rays (center and ±1/3 width perpendicular). Filter hits by minimum Z gap to exclude same-layer intersections. Interpolate Z along hit segment (fixes tilted quad geometry). Average valid gaps to compute multiplier = avgGap / nomGap.

**Key Data Flow:**
- Input: `DataTree<Line>` segments preserving any branch structure (spirals, layers, mixed)
- Processing: Flattened for parallel computation, branch info stored for reconstruction
- Output: `DataTree<double>` multipliers matching input tree structure

## Code Conventions

**Grasshopper Patterns:**
- Use `DataTree<T>` for branched data, reconstruct via `GH_Path`
- `GH_ScriptInstance` is the base class for script components

**Performance Patterns:**
- `Parallel.For` for independent computations
- Pre-allocate arrays before parallel loops
- Use `lock` for aggregate statistics, `Interlocked` for counters
- `Mesh.Compact()` after construction

**Numerical Precision:**
- Use `1e-10` threshold for degenerate geometry checks
- Clamp interpolation parameters to [0,1]
- Handle vertical/degenerate segments explicitly

**Dependencies (referenced via Grasshopper, not NuGet):**
- `Rhino.Geometry` - 3D primitives and operations
- `Rhino.Geometry.Intersect` - Ray-mesh intersection
- `Grasshopper.Kernel.Data` - DataTree structures
