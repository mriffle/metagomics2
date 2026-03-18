/**
 * Centralized color definitions for light and dark mode.
 *
 * All visualization colors (taxonomy palettes, GO DAG node/edge colors,
 * status badge colors, etc.) are defined here so they can be imported by
 * any component that needs them.  CSS custom properties in index.css
 * mirror the same values for elements styled purely through CSS.
 */

// ---------------------------------------------------------------------------
// Taxonomy domain palettes
// ---------------------------------------------------------------------------

/** Light-mode domain palette — colorblind-friendly, medium saturation */
export const DOMAIN_PALETTE_LIGHT: number[][] = [
  [59, 130, 246],   // blue-500
  [16, 185, 129],   // emerald-500
  [245, 158, 11],   // amber-500
  [239, 68, 68],    // red-500
  [139, 92, 246],   // violet-500
  [236, 72, 153],   // pink-500
  [14, 165, 233],   // sky-500
  [234, 179, 8],    // yellow-500
  [168, 85, 247],   // purple-500
  [20, 184, 166],   // teal-500
]

/** Dark-mode domain palette — vivid, high-saturation neon */
export const DOMAIN_PALETTE_DARK: number[][] = [
  [96, 165, 250],   // bright blue
  [52, 211, 153],   // bright emerald
  [251, 191, 36],   // bright amber
  [248, 113, 113],  // bright red
  [167, 139, 250],  // bright violet
  [244, 114, 182],  // bright pink
  [56, 189, 248],   // bright sky
  [250, 204, 21],   // bright yellow
  [192, 132, 252],  // bright purple
  [45, 212, 191],   // bright teal
]

// ---------------------------------------------------------------------------
// Taxonomy chart colors
// ---------------------------------------------------------------------------

export const TAX_CHART = {
  light: {
    rootColor: 'rgb(255,255,255)',
    fallbackColor: 'rgb(229,231,235)',
    lineColor: '#ffffff',
    textColor: undefined as string | undefined,
    pathbarTextColor: undefined as string | undefined,
  },
  dark: {
    rootColor: 'rgb(17,24,39)',
    fallbackColor: 'rgb(55,65,81)',
    lineColor: '#1f2937',
    textColor: '#111827',
    pathbarTextColor: '#e5e7eb',
  },
} as const

// ---------------------------------------------------------------------------
// GO DAG Cytoscape colors
// ---------------------------------------------------------------------------

export const GO_DAG = {
  light: {
    nodeBorderWidth: 1,
    nodeBorderColor: '#94a3b8',
    nodeBgColor: '#e2e8f0',
    nodeTextColor: '#1e1b4b',
    nodeShadow: null as null | { blur: number; color: string },

    edgeColor: '#cbd5e1',

    highlightBorderWidth: 2,
    highlightBorderColor: '#4338ca',
    highlightShadow: null as null | { blur: number; color: string },

    highlightEdgeColor: '#818cf8',

    exportBg: '#ffffff',

    colorRampDim: [255, 255, 255] as number[],
    colorRampDark: [0, 0, 0] as number[],
    colorRampDimBlend: 0.07,
    textLight: '#ffffff',
    textDark: '#1e1b4b',
  },
  dark: {
    nodeBorderWidth: 1.5,
    nodeBorderColor: '#6366f1',
    nodeBgColor: '#1e1b4b',
    nodeTextColor: '#e0e7ff',
    nodeShadow: { blur: 12, color: 'rgba(99,102,241,0.4)' },

    edgeColor: '#4b5563',

    highlightBorderWidth: 3,
    highlightBorderColor: '#a5b4fc',
    highlightShadow: { blur: 20, color: 'rgba(165,180,252,0.6)' },

    highlightEdgeColor: '#a5b4fc',

    exportBg: '#030712',

    colorRampDim: [30, 30, 40] as number[],
    colorRampDark: [255, 255, 255] as number[],
    colorRampDimBlend: 0.7,
    textLight: '#e0e7ff',
    textDark: '#1e1b4b',
  },
} as const

// ---------------------------------------------------------------------------
// Sankey-specific overrides
// ---------------------------------------------------------------------------

export const SANKEY = {
  light: {
    nodeLineColor: '#888',
  },
  dark: {
    nodeLineColor: '#4b5563',
  },
} as const

// ---------------------------------------------------------------------------
// Plotly layout colors
// ---------------------------------------------------------------------------

export const PLOTLY_LAYOUT = {
  light: {
    fontColor: undefined as string | undefined,
  },
  dark: {
    fontColor: '#e5e7eb',
  },
} as const

// ---------------------------------------------------------------------------
// Status badge classes (used by JobPage, AdminPage, HomePage)
// ---------------------------------------------------------------------------

export const STATUS_BADGE_CLASSES: Record<string, string> = {
  uploaded: 'bg-gray-100 text-gray-800 dark:bg-gray-500/20 dark:text-gray-300 dark:ring-1 dark:ring-gray-500/40',
  queued:   'bg-yellow-100 text-yellow-800 dark:bg-yellow-500/20 dark:text-yellow-300 dark:ring-1 dark:ring-yellow-500/40',
  pending:  'bg-yellow-100 text-yellow-800 dark:bg-yellow-500/20 dark:text-yellow-300 dark:ring-1 dark:ring-yellow-500/40',
  running:  'bg-blue-100 text-blue-800 dark:bg-blue-500/20 dark:text-blue-300 dark:ring-1 dark:ring-blue-500/40',
  completed:'bg-green-100 text-green-800 dark:bg-green-500/20 dark:text-green-300 dark:ring-1 dark:ring-green-500/40',
  done:     'bg-green-100 text-green-800 dark:bg-green-500/20 dark:text-green-300 dark:ring-1 dark:ring-green-500/40',
  failed:   'bg-red-100 text-red-800 dark:bg-red-500/20 dark:text-red-300 dark:ring-1 dark:ring-red-500/40',
}

/** Fallback badge class when the status key is not found */
export const STATUS_BADGE_DEFAULT =
  'bg-gray-100 text-gray-800 dark:bg-gray-500/20 dark:text-gray-300 dark:ring-1 dark:ring-gray-500/40'
