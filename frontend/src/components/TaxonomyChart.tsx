import { useCallback, useMemo, useRef, useState, useLayoutEffect } from 'react'
import createPlotlyComponent from 'react-plotly.js/factory'
import Plotly from 'plotly.js-dist-min'
import type { TaxonNode } from '../utils/taxonomyParser'
import type { ChartType } from '../pages/TaxonomyPage'

const Plot = createPlotlyComponent(Plotly)

// Categorical palette for domains — visually distinct, colorblind-friendly
const DOMAIN_PALETTE = [
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

/**
 * Blend a domain base color toward white based on depth ratio.
 * depth 0 (domain itself) = full color, depth 1 (species) = lightest.
 */
function domainColorAtDepth(base: number[], depthRatio: number): string {
  const t = Math.max(0, Math.min(1, depthRatio))
  // Blend from base color (t=0) toward a very light tint (t=1)
  const lightness = 0.35 + t * 0.55  // range [0.35, 0.90] — never fully white
  const r = Math.round(base[0] + (255 - base[0]) * lightness)
  const g = Math.round(base[1] + (255 - base[1]) * lightness)
  const b = Math.round(base[2] + (255 - base[2]) * lightness)
  return `rgb(${r},${g},${b})`
}

interface TaxonomyChartProps {
  nodes: TaxonNode[]
  chartType: ChartType
}

export default function TaxonomyChart({ nodes, chartType }: TaxonomyChartProps) {
  const containerRef = useRef<HTMLDivElement>(null)
  const tooltipRef = useRef<HTMLDivElement>(null)
  const [tooltip, setTooltip] = useState<{
    x: number
    y: number
    node: TaxonNode
  } | null>(null)
  const [tooltipPos, setTooltipPos] = useState<{ left: number; top: number }>({ left: 0, top: 0 })

  // Reposition tooltip to stay within viewport
  useLayoutEffect(() => {
    if (!tooltip || !tooltipRef.current) return
    const tRect = tooltipRef.current.getBoundingClientRect()

    let left = tooltip.x + 16
    let top = tooltip.y - 16

    if (left + tRect.width > window.innerWidth - 8) {
      left = tooltip.x - tRect.width - 16
    }
    if (left < 8) left = 8

    if (top + tRect.height > window.innerHeight - 8) {
      top = tooltip.y - tRect.height - 16
    }
    if (top < 8) top = 8

    setTooltipPos({ left, top })
  }, [tooltip])

  // Build a lookup map for tooltip data
  const nodeMap = useMemo(() => {
    const map = new Map<string, TaxonNode>()
    for (const n of nodes) map.set(n.taxId, n)
    return map
  }, [nodes])

  // Build Plotly data arrays
  const plotData = useMemo(() => {
    if (nodes.length === 0) return null

    // Find the root node(s) — nodes whose parentTaxId is empty or not in the set
    const nodeIds = new Set(nodes.map(n => n.taxId))

    const ids: string[] = []
    const labels: string[] = []
    const parents: string[] = []
    const values: number[] = []
    const colors: string[] = []
    const customdata: TaxonNode[] = []

    // Build parent lookup and find domain ancestor for each node
    const nodeById = new Map<string, TaxonNode>()
    for (const n of nodes) nodeById.set(n.taxId, n)

    // Walk up to find the domain ancestor for a given node
    function getDomainAncestor(node: TaxonNode): string | null {
      let current: TaxonNode | undefined = node
      const visited = new Set<string>()
      while (current) {
        if (visited.has(current.taxId)) break
        visited.add(current.taxId)
        if (current.rank === 'domain') return current.taxId
        current = current.parentTaxId ? nodeById.get(current.parentTaxId) : undefined
      }
      return null
    }

    // Collect distinct domains and assign colors
    const domainIds = new Set<string>()
    for (const n of nodes) {
      if (n.rank === 'domain') domainIds.add(n.taxId)
    }
    const domainList = Array.from(domainIds)
    const domainColorMap = new Map<string, number[]>()
    domainList.forEach((id, i) => {
      domainColorMap.set(id, DOMAIN_PALETTE[i % DOMAIN_PALETTE.length])
    })

    // Rank depth for lightness gradient
    const RANK_DEPTH: Record<string, number> = {
      domain: 0, kingdom: 1, phylum: 2, class: 3,
      order: 4, family: 5, genus: 6, species: 7,
    }
    const maxDepth = 7

    for (const node of nodes) {
      ids.push(node.taxId)
      labels.push(node.name)
      parents.push(nodeIds.has(node.parentTaxId) ? node.parentTaxId : '')
      values.push(node.quantity)
      customdata.push(node)

      if (node.rank === 'root') {
        colors.push('rgb(255,255,255)')
      } else {
        const domainId = getDomainAncestor(node)
        const base = domainId ? domainColorMap.get(domainId) : undefined
        if (base) {
          const depth = RANK_DEPTH[node.rank] ?? 0
          const depthRatio = depth / maxDepth
          colors.push(domainColorAtDepth(base, depthRatio))
        } else {
          colors.push('rgb(229,231,235)')  // gray-200 fallback
        }
      }
    }

    return { ids, labels, parents, values, colors, customdata }
  }, [nodes])

  // Keep nodeMap in a ref so hover callbacks don't need it as a dependency
  const nodeMapRef = useRef(nodeMap)
  nodeMapRef.current = nodeMap

  // Memoize Plotly props so tooltip state changes don't reset the chart
  const plotlyData = useMemo(() => {
    if (!plotData) return null
    const trace: any = {
      type: chartType,
      ids: plotData.ids,
      labels: plotData.labels,
      parents: plotData.parents,
      values: plotData.values,
      marker: {
        colors: plotData.colors,
        line: { width: chartType === 'treemap' ? 1 : 0.5, color: '#ffffff' },
      },
      branchvalues: 'total',
      hoverinfo: 'none',
      textinfo: 'label',
      textfont: { size: chartType === 'sunburst' ? 11 : 12 },
      insidetextorientation: 'auto',
      maxdepth: -1,
    }
    if (chartType === 'treemap') {
      trace.tiling = { packing: 'squarify', pad: 2 }
      trace.pathbar = { visible: true, textfont: { size: 12 } }
    }
    return [trace]
  }, [plotData, chartType])

  const plotlyLayout = useMemo(() => {
    const layout: any = {
      margin: { t: chartType === 'treemap' ? 20 : 10, b: 10, l: 10, r: 10 },
      paper_bgcolor: 'transparent',
      font: { family: 'Inter, system-ui, sans-serif' },
      autosize: true,
    }
    if (chartType === 'sunburst') {
      layout.sunburstcolorway = []
    }
    return layout
  }, [chartType])

  const plotlyConfig = useMemo(() => ({
    responsive: true,
    displayModeBar: false,
    displaylogo: false,
    modeBarButtonsToRemove: [
      'zoom2d', 'pan2d', 'select2d', 'lasso2d', 'zoomIn2d', 'zoomOut2d',
      'autoScale2d', 'hoverClosestCartesian', 'hoverCompareCartesian',
      'toggleSpikelines',
    ] as any[],
  }), [])

  const handleHover = useCallback((event: any) => {
    if (event.points && event.points.length > 0) {
      const point = event.points[0]
      const taxId = point.id
      const node = nodeMapRef.current.get(taxId)
      if (node && event.event) {
        setTooltip({
          x: event.event.clientX,
          y: event.event.clientY,
          node,
        })
      }
    }
  }, [])

  const handleUnhover = useCallback(() => {
    setTooltip(null)
  }, [])

  if (!plotlyData || nodes.length === 0) {
    return (
      <div className="flex items-center justify-center h-full text-gray-500">
        No taxonomy data available.
      </div>
    )
  }

  return (
    <div ref={containerRef} className="relative w-full h-full">
      <Plot
        data={plotlyData}
        layout={plotlyLayout}
        config={plotlyConfig}
        style={{ width: '100%', height: '100%' }}
        useResizeHandler
        onHover={handleHover}
        onUnhover={handleUnhover}
      />

      {/* Tooltip */}
      {tooltip && (
        <div
          ref={tooltipRef}
          className="fixed z-50 bg-white border border-gray-300 rounded-lg shadow-lg p-3 text-sm pointer-events-none"
          style={{
            left: tooltipPos.left,
            top: tooltipPos.top,
            maxWidth: 320,
          }}
        >
          <p className="font-mono text-xs text-gray-500 mb-1">Tax ID: {tooltip.node.taxId}</p>
          <p className="font-semibold text-gray-900 mb-1">{tooltip.node.name}</p>
          <p className="text-xs text-gray-500 mb-2 capitalize">{tooltip.node.rank}</p>
          <div className="grid grid-cols-2 gap-x-4 gap-y-1 text-xs">
            <span className="text-gray-600">Quantity:</span>
            <span className="font-medium text-right">{tooltip.node.quantity > 9999 ? tooltip.node.quantity.toExponential(3) : tooltip.node.quantity.toFixed(2)}</span>
            <span className="text-gray-600">Ratio (Total):</span>
            <span className="font-medium text-right">{(tooltip.node.ratioTotal * 100).toFixed(4)}%</span>
            <span className="text-gray-600">Ratio (Annotated):</span>
            <span className="font-medium text-right">{(tooltip.node.ratioAnnotated * 100).toFixed(4)}%</span>
            <span className="text-gray-600"># Peptides:</span>
            <span className="font-medium text-right">{tooltip.node.nPeptides}</span>
          </div>
        </div>
      )}
    </div>
  )
}
