//
//  MoleculeDrawing.swift
//  MacQC
//
//  ΨLab 5.1: drawings of the molecule. The Optimisation tab draws the geometry of every
//  optimisation step (bond lengths and the H–O–H angle written on it, the starting geometry
//  dashed) and plays the steps one after the other; the Frequencies tab animates every normal
//  mode. H₂, LiH and H₂O are linear or planar, so a drawing in the plane in which the atoms
//  spread shows every atom and every displacement without distortion.
//

import SwiftUI
import AppKit

// MARK: - Geometry helpers

/// Cartesian coordinates (bohr) of the atoms, three per atom.
func atomCoordinates(_ atoms: [MacQCAtom]) -> [[Double]] {
    atoms.map { a in (0..<3).map { $0 < a.position.count ? a.position[$0] : 0 } }
}

/// Index of the heaviest atom (the first of equal ones).
func heaviestAtom(_ atoms: [MacQCAtom]) -> Int? {
    guard !atoms.isEmpty else { return nil }
    var hub = 0
    for (i, a) in atoms.enumerated() where a.Z > atoms[hub].Z { hub = i }
    return hub
}

/// Bonds of the molecules of ΨLab: every atom bonded to the heaviest one (H₂: H–H; LiH: Li–H;
/// H₂O: the two O–H bonds and no H–H bond), whatever the bond lengths during an optimisation.
func moleculeBonds(_ atoms: [MacQCAtom]) -> [(Int, Int)] {
    guard atoms.count > 1, let hub = heaviestAtom(atoms) else { return [] }
    return atoms.indices.filter { $0 != hub }.map { (min(hub, $0), max(hub, $0)) }
}

func vectorDistance(_ a: [Double], _ b: [Double]) -> Double {
    var s = 0.0
    for k in 0..<min(a.count, b.count) { s += (a[k] - b[k]) * (a[k] - b[k]) }
    return s.squareRoot()
}

func vectorNorm(_ a: [Double]) -> Double {
    a.reduce(0) { $0 + $1 * $1 }.squareRoot()
}

/// Positions moved along a displacement: r + factor · d (missing displacements count as zero).
func displacedCoordinates(_ positions: [[Double]], _ displacement: [[Double]], _ factor: Double) -> [[Double]] {
    positions.enumerated().map { i, p in
        let d = i < displacement.count && displacement[i].count == 3 ? displacement[i] : [0, 0, 0]
        return (0..<3).map { p[$0] + factor * d[$0] }
    }
}

// MARK: - Projection and drawing

/// The view on a set of geometries of one molecule: the two Cartesian axes along which the atoms
/// spread most (horizontal and vertical), the centre and the spans of the region to show. The
/// heaviest atom is put on top (H₂O is drawn as the usual inverted V).
struct MoleculeFrame {
    var horizontal = 2
    var vertical = 0
    var flipVertical = false
    var centerH = 0.0
    var centerV = 0.0
    var spreadH = 1.0
    var spreadV = 0.0

    init(geometries: [[[Double]]], heavyAtom: Int?) {
        let all = geometries.flatMap { $0 }
        guard !all.isEmpty else { return }
        var lo = [Double](repeating: Double.infinity, count: 3)
        var hi = [Double](repeating: -Double.infinity, count: 3)
        for p in all {
            for k in 0..<3 {
                lo[k] = min(lo[k], p[k])
                hi[k] = max(hi[k], p[k])
            }
        }
        let spread = (0..<3).map { hi[$0] - lo[$0] }
        let order = [0, 1, 2].sorted { spread[$0] > spread[$1] + 1e-9 || (abs(spread[$0] - spread[$1]) <= 1e-9 && $0 < $1) }
        horizontal = order[0]
        vertical = order[1]
        centerH = (hi[horizontal] + lo[horizontal]) / 2
        centerV = (hi[vertical] + lo[vertical]) / 2
        spreadH = spread[horizontal]
        spreadV = spread[vertical]
        if let h = heavyAtom, let first = geometries.first, h < first.count, first.count > 1 {
            let others = first.indices.filter { $0 != h }.map { first[$0][vertical] }
            let mean = others.reduce(0, +) / Double(others.count)
            flipVertical = first[h][vertical] < mean - 1e-6
        }
    }

    /// Pixels per bohr so that the region fits into `size` with `padding` pixels around it.
    func scale(for size: CGSize, padding: Double) -> Double {
        let w = max(Double(size.width) - 2 * padding, 10)
        let h = max(Double(size.height) - 2 * padding, 10)
        var s = w / max(spreadH, 0.5)
        if spreadV > 1e-6 { s = min(s, h / spreadV) }
        return min(s, 150)
    }

    func point(_ p: [Double], size: CGSize, scale: Double) -> CGPoint {
        let x = Double(size.width) / 2 + (p[horizontal] - centerH) * scale
        let dy = (p[vertical] - centerV) * scale
        return CGPoint(x: x, y: Double(size.height) / 2 + (flipVertical ? dy : -dy))
    }
}

/// Colours and sizes of the atoms (CPK colours: H white, He cyan, Li violet, O red).
enum AtomStyle {
    static func color(_ Z: Int) -> Color {
        switch Z {
        case 1: return Color(white: 0.93)
        case 2: return Color(red: 0.75, green: 0.95, blue: 1.0)
        case 3: return Color(red: 0.75, green: 0.45, blue: 1.0)
        case 8: return Color(red: 0.92, green: 0.16, blue: 0.12)
        default: return Color.gray
        }
    }

    static func radius(_ Z: Int) -> Double {
        switch Z {
        case 1, 2: return 11
        case 3: return 17
        default: return 15
        }
    }

    static func labelColor(_ Z: Int) -> Color {
        Z <= 2 ? Color.black : Color.white
    }
}

/// Draws atoms, bonds, outlines, arrows and labels into a Canvas for one MoleculeFrame.
struct MoleculeSketch {
    let frame: MoleculeFrame
    let size: CGSize
    let scale: Double

    init(frame: MoleculeFrame, size: CGSize, padding: Double) {
        self.frame = frame
        self.size = size
        self.scale = frame.scale(for: size, padding: padding)
    }

    func point(_ p: [Double]) -> CGPoint {
        frame.point(p, size: size, scale: scale)
    }

    func bonds(_ gc: GraphicsContext, _ positions: [[Double]], _ bonds: [(Int, Int)], color: Color,
               width: Double, dashed: Bool = false) {
        for (i, j) in bonds where i < positions.count && j < positions.count {
            var path = Path()
            path.move(to: point(positions[i]))
            path.addLine(to: point(positions[j]))
            gc.stroke(path, with: .color(color),
                      style: StrokeStyle(lineWidth: width, lineCap: .round, dash: dashed ? [5, 4] : []))
        }
    }

    func atoms(_ gc: GraphicsContext, _ atoms: [MacQCAtom], _ positions: [[Double]]) {
        for (i, atom) in atoms.enumerated() where i < positions.count {
            let p = point(positions[i])
            let r = AtomStyle.radius(atom.Z)
            let rect = CGRect(x: Double(p.x) - r, y: Double(p.y) - r, width: 2 * r, height: 2 * r)
            let highlight = CGPoint(x: Double(p.x) - 0.35 * r, y: Double(p.y) - 0.35 * r)
            gc.fill(Path(ellipseIn: rect),
                    with: .radialGradient(Gradient(colors: [Color.white, AtomStyle.color(atom.Z)]),
                                          center: highlight, startRadius: 0, endRadius: 1.3 * r))
            gc.stroke(Path(ellipseIn: rect), with: .color(Color.primary.opacity(0.55)), lineWidth: 1)
            gc.draw(Text(atom.symbol).font(.caption.bold()).foregroundStyle(AtomStyle.labelColor(atom.Z)), at: p)
        }
    }

    /// Dashed circles: a reference geometry (equilibrium, or the start of an optimisation).
    func outlines(_ gc: GraphicsContext, _ atoms: [MacQCAtom], _ positions: [[Double]]) {
        for (i, atom) in atoms.enumerated() where i < positions.count {
            let p = point(positions[i])
            let r = AtomStyle.radius(atom.Z)
            let rect = CGRect(x: Double(p.x) - r, y: Double(p.y) - r, width: 2 * r, height: 2 * r)
            gc.stroke(Path(ellipseIn: rect), with: .color(Color.secondary.opacity(0.8)),
                      style: StrokeStyle(lineWidth: 1.2, dash: [4, 3]))
        }
    }

    func arrow(_ gc: GraphicsContext, from a: [Double], to b: [Double], color: Color) {
        let base = point(a)
        let tip = point(b)
        let dx = Double(tip.x - base.x)
        let dy = Double(tip.y - base.y)
        let length = (dx * dx + dy * dy).squareRoot()
        guard length > 6 else { return }
        let ux = dx / length
        let uy = dy / length
        let tx = Double(tip.x)
        let ty = Double(tip.y)
        var path = Path()
        path.move(to: base)
        path.addLine(to: tip)
        path.move(to: CGPoint(x: tx - 9 * ux + 5 * uy, y: ty - 9 * uy - 5 * ux))
        path.addLine(to: tip)
        path.addLine(to: CGPoint(x: tx - 9 * ux - 5 * uy, y: ty - 9 * uy + 5 * ux))
        gc.stroke(path, with: .color(color), style: StrokeStyle(lineWidth: 2.2, lineCap: .round, lineJoin: .round))
    }

    /// A small text on a background, centred at `p`.
    func label(_ gc: GraphicsContext, _ string: String, at p: CGPoint) {
        let text = gc.resolve(Text(string).font(.caption.monospacedDigit()).foregroundStyle(Color.primary))
        let s = text.measure(in: CGSize(width: 300, height: 40))
        let rect = CGRect(x: Double(p.x) - Double(s.width) / 2 - 4, y: Double(p.y) - Double(s.height) / 2 - 2,
                          width: Double(s.width) + 8, height: Double(s.height) + 4)
        gc.fill(Path(roundedRect: rect, cornerRadius: 4),
                with: .color(Color(nsColor: .controlBackgroundColor).opacity(0.85)))
        gc.draw(text, at: p, anchor: .center)
    }
}

// MARK: - Optimisation step

/// The molecule at one optimisation step, drawn to the scale of the whole optimisation (playing
/// the steps shows how it moves), with its bond lengths (Å) and, for H₂O, the H–O–H angle;
/// optionally the starting geometry as dashed outlines.
struct OptimizationGeometryCanvas: View {
    var atoms: [MacQCAtom]
    var start: [MacQCAtom]?
    var frame: MoleculeFrame

    var body: some View {
        Canvas { gc, size in
            draw(gc, size: size)
        }
        .accessibilityLabel("Molecule at this optimisation step")
    }

    private func draw(_ gc: GraphicsContext, size: CGSize) {
        guard !atoms.isEmpty else { return }
        let sketch = MoleculeSketch(frame: frame, size: size, padding: 40)
        let positions = atomCoordinates(atoms)
        let bonds = moleculeBonds(atoms)
        if let start, start.count == atoms.count {
            let startPositions = atomCoordinates(start)
            sketch.bonds(gc, startPositions, moleculeBonds(start), color: Color.secondary.opacity(0.7), width: 1.5, dashed: true)
            sketch.outlines(gc, start, startPositions)
        }
        sketch.bonds(gc, positions, bonds, color: Color.secondary, width: 5)
        if bonds.count == 2, let hub = heaviestAtom(atoms) {
            angleMark(gc, sketch, positions, hub: hub, ends: bonds.map { $0.0 == hub ? $0.1 : $0.0 })
        }
        sketch.atoms(gc, atoms, positions)
        // bond lengths beside each bond, on the side away from the centre of the molecule
        let points = positions.map { sketch.point($0) }
        let cx = points.map { Double($0.x) }.reduce(0, +) / Double(points.count)
        let cy = points.map { Double($0.y) }.reduce(0, +) / Double(points.count)
        for (i, j) in bonds {
            let a = points[i], b = points[j]
            let dx = Double(b.x - a.x), dy = Double(b.y - a.y)
            let length = max((dx * dx + dy * dy).squareRoot(), 1)
            var nx = -dy / length, ny = dx / length
            let mx = (Double(a.x) + Double(b.x)) / 2, my = (Double(a.y) + Double(b.y)) / 2
            let side = nx * (cx - mx) + ny * (cy - my)
            if side > 1e-6 || (abs(side) <= 1e-6 && ny > 0) {
                nx = -nx
                ny = -ny
            }
            let r = vectorDistance(positions[i], positions[j]) * Units.bohrToAngstrom
            sketch.label(gc, String(format: "%.4f Å", r), at: CGPoint(x: mx + 18 * nx, y: my + 18 * ny))
        }
    }

    /// Arc and value of the angle at the central atom.
    private func angleMark(_ gc: GraphicsContext, _ sketch: MoleculeSketch, _ positions: [[Double]], hub: Int, ends: [Int]) {
        guard ends.count == 2 else { return }
        let u = (0..<3).map { positions[ends[0]][$0] - positions[hub][$0] }
        let v = (0..<3).map { positions[ends[1]][$0] - positions[hub][$0] }
        let cosine = (0..<3).map { u[$0] * v[$0] }.reduce(0, +) / max(vectorNorm(u) * vectorNorm(v), 1e-12)
        let degrees = acos(min(1, max(-1, cosine))) * 180 / Double.pi
        let c = sketch.point(positions[hub])
        let p1 = sketch.point(positions[ends[0]]), p2 = sketch.point(positions[ends[1]])
        let a1 = atan2(Double(p1.y - c.y), Double(p1.x - c.x))
        let a2 = atan2(Double(p2.y - c.y), Double(p2.x - c.x))
        var diff = a2 - a1
        while diff > Double.pi { diff -= 2 * Double.pi }
        while diff <= -Double.pi { diff += 2 * Double.pi }
        let radius = 30.0
        var arc = Path()
        for k in 0...24 {
            let t = a1 + diff * Double(k) / 24
            let q = CGPoint(x: Double(c.x) + radius * cos(t), y: Double(c.y) + radius * sin(t))
            if k == 0 { arc.move(to: q) } else { arc.addLine(to: q) }
        }
        gc.stroke(arc, with: .color(Color.accentColor), lineWidth: 1.5)
        let mid = a1 + diff / 2
        sketch.label(gc, String(format: "%.2f°", degrees),
                     at: CGPoint(x: Double(c.x) + 52 * cos(mid), y: Double(c.y) + 52 * sin(mid)))
    }
}

// MARK: - Normal mode

/// One normal mode as an animation: the atoms oscillate about the equilibrium geometry (dashed)
/// along the mode's Cartesian displacement, at a speed proportional to the wavenumber;
/// optionally with arrows to the turning points.
struct NormalModeCanvas: View {
    var atoms: [MacQCAtom]
    var displacement: [[Double]]
    var wavenumber: Double
    var animate: Bool
    /// 1: the largest atomic displacement is a quarter of the shortest bond.
    var amplitude: Double
    var arrows: Bool

    /// Seconds per vibration: 2 s at 1600 cm⁻¹, proportionally faster for higher wavenumbers.
    private var period: Double {
        min(2.4, max(0.6, 2.0 * 1600.0 / max(abs(wavenumber), 1.0)))
    }

    var body: some View {
        TimelineView(.animation(minimumInterval: 1.0 / 30.0, paused: !animate)) { context in
            let phase = animate ? sin(context.date.timeIntervalSinceReferenceDate * 2 * Double.pi / period) : 0
            Canvas { gc, size in
                draw(gc, size: size, phase: phase)
            }
        }
    }

    /// Bohr per unit of the normalised displacement.
    private func bohrPerUnit(_ positions: [[Double]]) -> Double {
        let shortest = moleculeBonds(atoms).map { vectorDistance(positions[$0.0], positions[$0.1]) }.min() ?? 1.8
        let largest = displacement.map { vectorNorm($0) }.max() ?? 1
        return amplitude * 0.25 * shortest / max(largest, 1e-6)
    }

    private func draw(_ gc: GraphicsContext, size: CGSize, phase: Double) {
        guard !atoms.isEmpty else { return }
        let equilibrium = atomCoordinates(atoms)
        let a = bohrPerUnit(equilibrium)
        let frame = MoleculeFrame(geometries: [displacedCoordinates(equilibrium, displacement, a),
                                               displacedCoordinates(equilibrium, displacement, -a)],
                                  heavyAtom: heaviestAtom(atoms))
        let sketch = MoleculeSketch(frame: frame, size: size, padding: 28)
        let bonds = moleculeBonds(atoms)
        sketch.bonds(gc, equilibrium, bonds, color: Color.secondary.opacity(0.5), width: 1.5, dashed: true)
        sketch.outlines(gc, atoms, equilibrium)
        if arrows {
            let turning = displacedCoordinates(equilibrium, displacement, a)
            for i in equilibrium.indices {
                sketch.arrow(gc, from: equilibrium[i], to: turning[i], color: Color.orange)
            }
        }
        let now = displacedCoordinates(equilibrium, displacement, a * phase)
        sketch.bonds(gc, now, bonds, color: Color.secondary, width: 5)
        sketch.atoms(gc, atoms, now)
    }
}
