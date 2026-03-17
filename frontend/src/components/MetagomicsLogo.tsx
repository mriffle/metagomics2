/**
 * MetagomicsLogo — Combined sunburst icon + stylized "Metagomics 2" wordmark.
 *
 * Single SVG that adapts to light/dark mode via CSS classes (Tailwind `dark:`
 * variant toggles the `.dark` class on `<html>`). The large "2" is rendered
 * behind "Metagomics" at reduced opacity, slanted counter-clockwise, and in a
 * contrasting accent color.
 */
export default function MetagomicsLogo({ className }: { className?: string }) {
  // Sunburst centred at (24, 24), r_outer = 22
  const CX = 24, CY = 24

  return (
    <svg
      xmlns="http://www.w3.org/2000/svg"
      viewBox="0 0 270 48"
      className={className}
      aria-label="Metagomics 2 logo"
      role="img"
    >
      <defs>
        <style>{`
          /* ---- Light mode (default) ---- */
          .mg-w-a { fill: #4f46e5; }
          .mg-w-b { fill: #818cf8; }
          .mg-w-c { fill: #059669; }
          .mg-w-d { fill: #34d399; }
          .mg-w-e { fill: #d97706; }
          .mg-w-f { fill: #fbbf24; }
          .mg-w-g { fill: #dc2626; }
          .mg-w-h { fill: #f87171; }
          .mg-w-i { fill: #7c3aed; }
          .mg-w-j { fill: #a78bfa; }
          .mg-center { fill: #312e81; }
          .mg-gap   { stroke: #ffffff; stroke-width: 0.8; }
          .mg-txt   { fill: #1e1b4b; }
          .mg-big2  { fill: #e11d48; }

          /* ---- Dark mode ---- */
          .dark .mg-w-a { fill: #818cf8; }
          .dark .mg-w-b { fill: #a5b4fc; }
          .dark .mg-w-c { fill: #34d399; }
          .dark .mg-w-d { fill: #6ee7b7; }
          .dark .mg-w-e { fill: #fbbf24; }
          .dark .mg-w-f { fill: #fcd34d; }
          .dark .mg-w-g { fill: #f87171; }
          .dark .mg-w-h { fill: #fca5a5; }
          .dark .mg-w-i { fill: #a78bfa; }
          .dark .mg-w-j { fill: #c4b5fd; }
          .dark .mg-center { fill: #a5b4fc; }
          .dark .mg-gap   { stroke: #111827; stroke-width: 0.8; }
          .dark .mg-txt   { fill: #e0e7ff; }
          .dark .mg-big2  { fill: #fb7185; }
        `}</style>
      </defs>

      {/* ==== Sunburst icon ==== */}
      <g className="mg-gap">
        {/* Outer ring — 10 wedges (36° each) */}
        <path className="mg-w-a" d={wedge(CX, CY, 14, 22, 0, 36)} />
        <path className="mg-w-b" d={wedge(CX, CY, 14, 22, 36, 72)} />
        <path className="mg-w-c" d={wedge(CX, CY, 14, 22, 72, 108)} />
        <path className="mg-w-d" d={wedge(CX, CY, 14, 22, 108, 144)} />
        <path className="mg-w-e" d={wedge(CX, CY, 14, 22, 144, 180)} />
        <path className="mg-w-f" d={wedge(CX, CY, 14, 22, 180, 216)} />
        <path className="mg-w-g" d={wedge(CX, CY, 14, 22, 216, 252)} />
        <path className="mg-w-h" d={wedge(CX, CY, 14, 22, 252, 288)} />
        <path className="mg-w-i" d={wedge(CX, CY, 14, 22, 288, 324)} />
        <path className="mg-w-j" d={wedge(CX, CY, 14, 22, 324, 360)} />

        {/* Inner ring — 5 wedges (72° each) */}
        <path className="mg-w-a" d={wedge(CX, CY, 7, 14, 0, 72)} />
        <path className="mg-w-c" d={wedge(CX, CY, 7, 14, 72, 144)} />
        <path className="mg-w-e" d={wedge(CX, CY, 7, 14, 144, 216)} />
        <path className="mg-w-g" d={wedge(CX, CY, 7, 14, 216, 288)} />
        <path className="mg-w-i" d={wedge(CX, CY, 7, 14, 288, 360)} />
      </g>
      {/* Centre dot (no gap stroke) */}
      <circle className="mg-center" cx={CX} cy={CY} r="6" />

      {/* ==== Wordmark ==== */}
      <g>
        {/* "Metagomics" — foreground, bold */}
        <text
          className="mg-txt"
          x="56"
          y="34"
          fontSize="26"
          fontWeight="800"
          fontFamily="'Inter','Helvetica Neue',Arial,sans-serif"
          letterSpacing="-0.5"
        >
          MetaGOmics
        </text>

        {/* Large "2" — punchy red, slanted CCW, after the word */}
        <text
          className="mg-big2"
          x="248"
          y="49"
          fontSize="62"
          fontWeight="900"
          fontFamily="'Inter','Helvetica Neue',Arial,sans-serif"
          transform="rotate(-12, 230, 26)"
          textAnchor="middle"
          dominantBaseline="auto"
        >
          2
        </text>
      </g>
    </svg>
  )
}

/* ---------- SVG arc math helper ---------- */

/** Build a `d` attribute for a single annular wedge (sector between two radii). */
function wedge(
  cx: number, cy: number,
  r1: number, r2: number,
  startDeg: number, endDeg: number,
): string {
  const rad = (d: number) => ((d - 90) * Math.PI) / 180
  const px = (r: number, d: number) => cx + r * Math.cos(rad(d))
  const py = (r: number, d: number) => cy + r * Math.sin(rad(d))
  const lg = endDeg - startDeg > 180 ? 1 : 0

  return [
    `M${px(r1, startDeg)},${py(r1, startDeg)}`,
    `L${px(r2, startDeg)},${py(r2, startDeg)}`,
    `A${r2},${r2} 0 ${lg} 1 ${px(r2, endDeg)},${py(r2, endDeg)}`,
    `L${px(r1, endDeg)},${py(r1, endDeg)}`,
    `A${r1},${r1} 0 ${lg} 0 ${px(r1, startDeg)},${py(r1, startDeg)}`,
    'Z',
  ].join(' ')
}
