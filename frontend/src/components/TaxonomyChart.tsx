import { useCallback, useMemo, useRef, useState, useLayoutEffect } from 'react'
import createPlotlyComponent from 'react-plotly.js/factory'
import Plotly from 'plotly.js-dist-min'
import type { TaxonNode } from '../utils/taxonomyParser'
import type { ChartType } from '../pages/TaxonomyPage'
import { useTheme } from '../ThemeContext'
import { DOMAIN_PALETTE_LIGHT, DOMAIN_PALETTE_DARK, TAX_CHART, SANKEY, PLOTLY_LAYOUT } from '../utils/colors'

const Plot = createPlotlyComponent(Plotly)


/**
 * Blend a domain base color toward white based on depth ratio.
 * depth 0 (domain itself) = full color, depth 1 (species) = lightest.
 */
function domainColorAtDepth(base: number[], depthRatio: number, isDark = false): string {
  const t = Math.max(0, Math.min(1, depthRatio))
  if (isDark) {
    // Dark mode: keep colors vivid and bright. Slight shift toward lighter tint at depth.
    // Range [0, 0.30] — minimal washout, stays saturated and neon.
    const shift = t * 0.30
    const r = Math.round(base[0] + (255 - base[0]) * shift)
    const g = Math.round(base[1] + (255 - base[1]) * shift)
    const b = Math.round(base[2] + (255 - base[2]) * shift)
    return `rgb(${r},${g},${b})`
  }
  // Light mode: blend from base color toward a very light tint
  const lightness = 0.35 + t * 0.55  // range [0.35, 0.90] — never fully white
  const r = Math.round(base[0] + (255 - base[0]) * lightness)
  const g = Math.round(base[1] + (255 - base[1]) * lightness)
  const b = Math.round(base[2] + (255 - base[2]) * lightness)
  return `rgb(${r},${g},${b})`
}

interface TaxonomyChartProps {
  nodes: TaxonNode[]
  chartType: ChartType
  filterLabel?: string
  onNodeClick?: (taxId: string | null) => void
}

export default function TaxonomyChart({ nodes, chartType, filterLabel, onNodeClick }: TaxonomyChartProps) {
  const containerRef = useRef<HTMLDivElement>(null)
  const tooltipRef = useRef<HTMLDivElement>(null)
  const { theme } = useTheme()
  const isDark = theme === 'dark'
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
    const ranks: string[] = []
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
    const palette = isDark ? DOMAIN_PALETTE_DARK : DOMAIN_PALETTE_LIGHT
    domainList.forEach((id, i) => {
      domainColorMap.set(id, palette[i % palette.length])
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
      ranks.push(node.rank)
      customdata.push(node)

      const tc = isDark ? TAX_CHART.dark : TAX_CHART.light
      if (node.rank === 'root') {
        colors.push(tc.rootColor)
      } else {
        const domainId = getDomainAncestor(node)
        const base = domainId ? domainColorMap.get(domainId) : undefined
        if (base) {
          const depth = RANK_DEPTH[node.rank] ?? 0
          const depthRatio = depth / maxDepth
          colors.push(domainColorAtDepth(base, depthRatio, isDark))
        } else {
          colors.push(tc.fallbackColor)
        }
      }
    }

    return { ids, labels, parents, values, colors, ranks, customdata }
  }, [nodes, isDark])

  // Keep nodeMap in a ref so hover callbacks don't need it as a dependency
  const nodeMapRef = useRef(nodeMap)
  nodeMapRef.current = nodeMap

  // Keep onNodeClick in a ref to avoid stale closures
  const onNodeClickRef = useRef(onNodeClick)
  onNodeClickRef.current = onNodeClick

  const handleClick = useCallback((event: any) => {
    if (chartType === 'sankey') return
    if (event.points && event.points.length > 0) {
      const point = event.points[0]
      const taxId = point.id
      if (taxId) onNodeClickRef.current?.(String(taxId))
    }
  }, [chartType])

  // Memoize Plotly props so tooltip state changes don't reset the chart
  const plotlyData = useMemo(() => {
    if (!plotData) return null

    // Sankey uses a completely different trace format (node + link arrays)
    if (chartType === 'sankey') {
      const idToIndex = new Map<string, number>()
      plotData.ids.forEach((id, i) => idToIndex.set(id, i))

      const source: number[] = []
      const target: number[] = []
      const value: number[] = []
      const linkColor: string[] = []

      for (let i = 0; i < plotData.ids.length; i++) {
        const parentId = plotData.parents[i]
        if (!parentId) continue
        const parentIdx = idToIndex.get(parentId)
        if (parentIdx === undefined) continue
        source.push(parentIdx)
        target.push(i)
        value.push(plotData.values[i])
        // Use a semi-transparent version of the target node color for links
        const c = plotData.colors[i]
        linkColor.push(c.replace('rgb(', 'rgba(').replace(')', ',0.4)'))
      }

      // Assign explicit x positions so each rank aligns to a column
      const RANK_COL: Record<string, number> = {
        root: 0, domain: 1, kingdom: 2, phylum: 3, class: 4,
        order: 5, family: 6, genus: 7, species: 8,
      }
      // Find the max column actually used to normalize x to [0, 1]
      let maxCol = 0
      for (const r of plotData.ranks) {
        const col = RANK_COL[r] ?? 0
        if (col > maxCol) maxCol = col
      }
      // Plotly requires x in (0, 1) exclusive — use small epsilon at edges
      const eps = 0.001
      const nodeX = plotData.ranks.map(r => {
        const col = RANK_COL[r] ?? 0
        return maxCol === 0 ? eps : eps + (col / maxCol) * (1 - 2 * eps)
      })

      return [{
        type: 'sankey',
        orientation: 'h',
        arrangement: 'snap',
        node: {
          label: plotData.labels,
          color: plotData.colors,
          x: nodeX,
          pad: 20,
          thickness: 15,
          line: { color: isDark ? SANKEY.dark.nodeLineColor : SANKEY.light.nodeLineColor, width: 0.5 },
        },
        link: { source, target, value, color: linkColor },
        valueformat: value.some(v => v > 9999) ? '.3e' : ',.2f',
        textfont: { color: isDark ? TAX_CHART.dark.textColor : TAX_CHART.light.textColor },
      }]
    }

    const trace: any = {
      type: chartType,
      ids: plotData.ids,
      labels: plotData.labels,
      parents: plotData.parents,
      values: plotData.values,
      marker: {
        colors: plotData.colors,
        line: { width: chartType === 'sunburst' ? 0.5 : 1, color: isDark ? TAX_CHART.dark.lineColor : TAX_CHART.light.lineColor },
      },
      branchvalues: 'total',
      hoverinfo: 'none',
      textinfo: 'label',
      textfont: { size: chartType === 'sunburst' ? 11 : 12, color: isDark ? TAX_CHART.dark.textColor : TAX_CHART.light.textColor },
      insidetextorientation: 'auto',
      maxdepth: -1,
    }
    if (chartType === 'treemap') {
      trace.tiling = { packing: 'squarify', pad: 2 }
      trace.pathbar = { visible: true, textfont: { size: 12, color: isDark ? TAX_CHART.dark.pathbarTextColor : TAX_CHART.light.pathbarTextColor } }
    }
    if (chartType === 'icicle') {
      trace.tiling = { orientation: 'v' }
      trace.pathbar = { visible: true, textfont: { size: 12, color: isDark ? TAX_CHART.dark.pathbarTextColor : TAX_CHART.light.pathbarTextColor } }
    }
    return [trace]
  }, [plotData, chartType, isDark])

  const plotlyLayout = useMemo(() => {
    const layout: any = {
      margin: { t: chartType === 'sunburst' ? 10 : 20, b: 10, l: 10, r: 10 },
      paper_bgcolor: 'transparent',
      plot_bgcolor: 'transparent',
      font: { family: 'Inter, system-ui, sans-serif', color: isDark ? PLOTLY_LAYOUT.dark.fontColor : PLOTLY_LAYOUT.light.fontColor },
      autosize: true,
    }
    if (chartType === 'sunburst') {
      layout.sunburstcolorway = []
    }
    return layout
  }, [chartType, isDark])

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
      <div className="flex items-center justify-center h-full text-gray-500 dark:text-gray-400">
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
        onHover={chartType !== 'sankey' ? handleHover : undefined}
        onUnhover={chartType !== 'sankey' ? handleUnhover : undefined}
        onClick={chartType !== 'sankey' ? handleClick : undefined}
      />

      {/* Tooltip (not used for Sankey — it has its own hover) */}
      {tooltip && chartType !== 'sankey' && (
        <div
          ref={tooltipRef}
          className="fixed z-50 bg-white dark:bg-gray-800 border border-gray-300 dark:border-gray-600 rounded-lg shadow-lg dark:shadow-indigo-500/10 p-3 text-sm pointer-events-none"
          style={{
            left: tooltipPos.left,
            top: tooltipPos.top,
            maxWidth: 320,
          }}
        >
          <p className="font-mono text-xs text-gray-500 dark:text-gray-400 mb-1">Tax ID: {tooltip.node.taxId}</p>
          <p className="font-semibold text-gray-900 dark:text-gray-100 mb-1">{tooltip.node.name}</p>
          <p className="text-xs text-gray-500 dark:text-gray-400 mb-2 capitalize">{tooltip.node.rank}</p>
          <div className="grid grid-cols-2 gap-x-4 gap-y-1 text-xs">
            <span className="text-gray-600 dark:text-gray-400">Quantity:</span>
            <span className="font-medium text-right text-gray-900 dark:text-gray-100">{tooltip.node.quantity > 9999 ? tooltip.node.quantity.toExponential(3) : tooltip.node.quantity.toFixed(2)}</span>
            <span className="text-gray-600 dark:text-gray-400">Ratio (Total):</span>
            <span className="font-medium text-right text-gray-900 dark:text-gray-100">{(tooltip.node.ratioTotal * 100).toFixed(4)}%</span>
            {filterLabel && tooltip.node.fractionOfTaxon != null && (
              <>
                <span className="text-gray-600 dark:text-gray-400">Fraction of Taxon:</span>
                <span className="font-medium text-right text-gray-900 dark:text-gray-100">{(tooltip.node.fractionOfTaxon * 100).toFixed(4)}%</span>
              </>
            )}
            {filterLabel && tooltip.node.fractionOfGo != null && (
              <>
                <span className="text-gray-600 dark:text-gray-400">Fraction of GO:</span>
                <span className="font-medium text-right text-gray-900 dark:text-gray-100">{(tooltip.node.fractionOfGo * 100).toFixed(4)}%</span>
              </>
            )}
            {!filterLabel && (
              <>
                <span className="text-gray-600 dark:text-gray-400">Ratio (Annotated):</span>
                <span className="font-medium text-right text-gray-900 dark:text-gray-100">{(tooltip.node.ratioAnnotated * 100).toFixed(4)}%</span>
              </>
            )}
            <span className="text-gray-600 dark:text-gray-400"># Peptides:</span>
            <span className="font-medium text-right text-gray-900 dark:text-gray-100">{tooltip.node.nPeptides}</span>
          </div>
        </div>
      )}
    </div>
  )
}
