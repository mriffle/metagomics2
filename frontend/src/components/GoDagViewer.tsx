import { useEffect, useLayoutEffect, useRef, useCallback, useState } from 'react'
import cytoscape from 'cytoscape'
import dagre from 'cytoscape-dagre'
import cytoscapeSvg from 'cytoscape-svg'
import type { GoTermNode, MetricKey } from '../pages/GoDagPage'
import { useTheme } from '../ThemeContext'
import { GO_DAG } from '../utils/colors'

// Register extensions
cytoscape.use(dagre)
cytoscape.use(cytoscapeSvg)

// Parse a hex color string to [r, g, b]
function hexToRgb(hex: string): number[] {
  const h = hex.replace('#', '')
  return [
    parseInt(h.substring(0, 2), 16),
    parseInt(h.substring(2, 4), 16),
    parseInt(h.substring(4, 6), 16),
  ]
}

// Blend two RGB colors by factor t (0 = a, 1 = b)
function blendRgb(a: number[], b: number[], t: number): number[] {
  return [
    Math.round(a[0] + (b[0] - a[0]) * t),
    Math.round(a[1] + (b[1] - a[1]) * t),
    Math.round(a[2] + (b[2] - a[2]) * t),
  ]
}

// Generate a 5-stop color ramp from near-white through baseColor to a darkened version
function generateColorStops(baseHex: string, isDark = false): number[][] {
  const base = hexToRgb(baseHex)
  const c = isDark ? GO_DAG.dark : GO_DAG.light
  if (isDark) {
    const dim = blendRgb(base, c.colorRampDim, c.colorRampDimBlend)
    const bright = blendRgb(base, c.colorRampDark, 0.3)
    return [
      dim,
      blendRgb(dim, base, 0.4),
      base,
      blendRgb(base, bright, 0.5),
      bright,
    ]
  }
  const white = c.colorRampDim
  const dark = blendRgb(base, c.colorRampDark, 0.6)
  return [
    blendRgb(white, base, c.colorRampDimBlend),
    blendRgb(white, base, 0.25),
    blendRgb(white, base, 0.55),
    base,
    dark,
  ]
}

function interpolateColor(t: number, stops: number[][]): string {
  const clamped = Math.max(0, Math.min(1, t))
  const segment = clamped * (stops.length - 1)
  const i = Math.floor(segment)
  const f = segment - i

  if (i >= stops.length - 1) {
    const c = stops[stops.length - 1]
    return `rgb(${c[0]},${c[1]},${c[2]})`
  }

  const c0 = stops[i]
  const c1 = stops[i + 1]
  const r = Math.round(c0[0] + (c1[0] - c0[0]) * f)
  const g = Math.round(c0[1] + (c1[1] - c0[1]) * f)
  const b = Math.round(c0[2] + (c1[2] - c0[2]) * f)
  return `rgb(${r},${g},${b})`
}

/**
 * Choose a readable text color based on the actual background color's luminance.
 * Uses WCAG relative-luminance formula instead of the abstract ramp position.
 */
function textColorForBg(bgColorStr: string, isDark = false): string {
  const c = isDark ? GO_DAG.dark : GO_DAG.light
  // Parse "rgb(r,g,b)" produced by interpolateColor
  const m = bgColorStr.match(/rgb\((\d+),(\d+),(\d+)\)/)
  if (!m) return c.textLight
  const [r, g, b] = [parseInt(m[1]), parseInt(m[2]), parseInt(m[3])]
  // sRGB → linear channel
  const toLinear = (v: number) => {
    const s = v / 255
    return s <= 0.03928 ? s / 12.92 : Math.pow((s + 0.055) / 1.055, 2.4)
  }
  const luminance = 0.2126 * toLinear(r) + 0.7152 * toLinear(g) + 0.0722 * toLinear(b)
  // Threshold ~0.35 works well: light backgrounds get dark text, dark backgrounds get light text
  return luminance > 0.35 ? c.textDark : c.textLight
}

function getMetricValue(node: GoTermNode, metric: MetricKey): number {
  switch (metric) {
    case 'quantity': return node.quantity
    case 'ratioTotal': return node.ratioTotal
    case 'ratioAnnotated': return node.ratioAnnotated
    case 'nPeptides': return node.nPeptides
    case 'fractionOfTaxon': return node.fractionOfTaxon ?? 0
    case 'fractionOfGo': return node.fractionOfGo ?? 0
    case 'qvalueGoForTaxon': return node.qvalueGoForTaxon ?? 1
  }
}

function normalizeValues(nodes: GoTermNode[], metric: MetricKey): Map<string, number> {
  const useLog = metric === 'quantity' || metric === 'nPeptides'
  const isQvalue = metric === 'qvalueGoForTaxon'
  const values = nodes.map(n => getMetricValue(n, metric))

  let transformed: number[]
  if (isQvalue) {
    // Invert: low q-value = high color intensity via -log10(q + eps)
    const eps = 1e-10
    transformed = values.map(v => -Math.log10(Math.max(v, eps)))
  } else if (useLog) {
    transformed = values.map(v => Math.log1p(v))
  } else {
    transformed = values
  }

  const min = Math.min(...transformed)
  const max = Math.max(...transformed)
  const range = max - min

  const result = new Map<string, number>()
  nodes.forEach((node, i) => {
    result.set(node.id, range > 0 ? (transformed[i] - min) / range : 0)
  })
  return result
}

interface GoDagViewerProps {
  nodes: GoTermNode[]
  metric: MetricKey
  filterLabel?: string
  baseColor?: string
  onNodeClick?: (nodeId: string | null) => void
}

export default function GoDagViewer({ nodes, metric, filterLabel, baseColor = '#4338ca', onNodeClick }: GoDagViewerProps) {
  const containerRef = useRef<HTMLDivElement>(null)
  const cyRef = useRef<cytoscape.Core | null>(null)
  const tooltipRef = useRef<HTMLDivElement>(null)
  const onNodeClickRef = useRef(onNodeClick)
  onNodeClickRef.current = onNodeClick
  const { theme } = useTheme()
  const isDark = theme === 'dark'
  const [tooltip, setTooltip] = useState<{
    x: number
    y: number
    node: GoTermNode
  } | null>(null)
  const [tooltipPos, setTooltipPos] = useState<{ left: number; top: number }>({ left: 0, top: 0 })

  // Reposition tooltip to stay within viewport
  useLayoutEffect(() => {
    if (!tooltip || !tooltipRef.current || !containerRef.current) return
    const wrapper = containerRef.current.parentElement!
    const wRect = wrapper.getBoundingClientRect()
    const tRect = tooltipRef.current.getBoundingClientRect()

    // Convert rendered position to viewport coordinates
    let left = wRect.left + tooltip.x + 12
    let top = wRect.top + tooltip.y - 12

    // Flip horizontally if overflowing right edge of viewport
    if (left + tRect.width > window.innerWidth - 8) {
      left = wRect.left + tooltip.x - tRect.width - 12
    }
    if (left < 8) left = 8

    // Flip vertically if overflowing bottom edge of viewport
    if (top + tRect.height > window.innerHeight - 8) {
      top = wRect.top + tooltip.y - tRect.height - 12
    }
    if (top < 8) top = 8

    setTooltipPos({ left, top })
  }, [tooltip])

  // Build cytoscape elements from nodes
  const buildElements = useCallback(() => {
    const nodeIds = new Set(nodes.map(n => n.id))
    const elements: cytoscape.ElementDefinition[] = []

    // Add nodes
    for (const node of nodes) {
      elements.push({
        group: 'nodes',
        data: {
          id: node.id,
          label: node.name.length > 40 ? node.name.substring(0, 37) + '...' : node.name,
          fullName: node.name,
          namespace: node.namespace,
          quantity: node.quantity,
          ratioTotal: node.ratioTotal,
          ratioAnnotated: node.ratioAnnotated,
          nPeptides: node.nPeptides,
        },
      })
    }

    // Add edges (parent -> child)
    for (const node of nodes) {
      for (const parentId of node.parentIds) {
        if (nodeIds.has(parentId)) {
          elements.push({
            group: 'edges',
            data: {
              id: `${parentId}->${node.id}`,
              source: parentId,
              target: node.id,
            },
          })
        }
      }
    }

    return elements
  }, [nodes])

  // Initialize cytoscape
  useEffect(() => {
    if (!containerRef.current || nodes.length === 0) return

    const elements = buildElements()
    const normalized = normalizeValues(nodes, metric)
    const stops = generateColorStops(baseColor, isDark)
    const colors = isDark ? GO_DAG.dark : GO_DAG.light

    // Build a lookup for quick access
    const nodeMap = new Map<string, GoTermNode>()
    for (const n of nodes) nodeMap.set(n.id, n)

    const cy = cytoscape({
      container: containerRef.current,
      elements,
      style: [
        {
          selector: 'node',
          style: {
            'shape': 'round-rectangle',
            'width': 'label',
            'height': 'label',
            'padding': '8px',
            'label': 'data(label)',
            'text-valign': 'center',
            'text-halign': 'center',
            'font-size': '12px',
            'text-wrap': 'wrap',
            'text-max-width': '140px',
            'border-width': colors.nodeBorderWidth,
            'border-color': colors.nodeBorderColor,
            'background-color': colors.nodeBgColor,
            'color': colors.nodeTextColor,
            'min-zoomed-font-size': 4,
            ...(colors.nodeShadow ? { 'shadow-blur': colors.nodeShadow.blur, 'shadow-color': colors.nodeShadow.color, 'shadow-offset-x': 0, 'shadow-offset-y': 0, 'shadow-opacity': 1 } : {}),
          } as any,
        },
        {
          selector: 'edge',
          style: {
            'width': 1,
            'line-color': colors.edgeColor,
            'target-arrow-color': colors.edgeColor,
            'target-arrow-shape': 'triangle',
            'arrow-scale': 0.6,
            'curve-style': 'bezier',
          } as any,
        },
        {
          selector: 'node:active',
          style: {
            'overlay-opacity': 0,
          } as any,
        },
        {
          selector: 'node.highlighted',
          style: {
            'border-width': colors.highlightBorderWidth,
            'border-color': colors.highlightBorderColor,
            ...(colors.highlightShadow ? { 'shadow-blur': colors.highlightShadow.blur, 'shadow-color': colors.highlightShadow.color, 'shadow-offset-x': 0, 'shadow-offset-y': 0, 'shadow-opacity': 1 } : {}),
          } as any,
        },
        {
          selector: 'edge.highlighted',
          style: {
            'width': 2,
            'line-color': colors.highlightEdgeColor,
            'target-arrow-color': colors.highlightEdgeColor,
          } as any,
        },
      ],
      layout: {
        name: 'dagre',
        rankDir: 'TB',
        ranker: 'tight-tree',
        nodeSep: 20,
        rankSep: 40,
        animate: false,
      } as any,
      minZoom: 0.05,
      maxZoom: 3,
      wheelSensitivity: 0.3,
    })

    // Apply colors based on metric
    cy.nodes().forEach((ele) => {
      const t = normalized.get(ele.id()) ?? 0
      const bgColor = interpolateColor(t, stops)
      ele.style('background-color', bgColor)
      ele.style('color', textColorForBg(bgColor, isDark))
    })

    // Collect all ancestor nodes + edges from a starting node up to the root(s)
    function collectAncestors(startNode: cytoscape.NodeSingular): { nodes: cytoscape.NodeCollection; edges: cytoscape.EdgeCollection } {
      let ancestorNodes = cy.collection() as cytoscape.NodeCollection
      let ancestorEdges = cy.collection() as cytoscape.EdgeCollection
      const visited = new Set<string>()
      const queue = [startNode]
      while (queue.length > 0) {
        const current = queue.shift()!
        if (visited.has(current.id())) continue
        visited.add(current.id())
        ancestorNodes = ancestorNodes.union(current) as cytoscape.NodeCollection
        // Incoming edges come from parent nodes (source=parent, target=current)
        const incoming = current.incomers('edge')
        ancestorEdges = ancestorEdges.union(incoming) as cytoscape.EdgeCollection
        incoming.forEach((edge: cytoscape.EdgeSingular) => {
          const src = edge.source()
          if (!visited.has(src.id())) queue.push(src)
        })
      }
      return { nodes: ancestorNodes, edges: ancestorEdges }
    }

    // Hover: show tooltip and highlight entire ancestor path to root
    cy.on('mouseover', 'node', (evt) => {
      const node = evt.target
      const { nodes: ancestorNodes, edges: ancestorEdges } = collectAncestors(node)
      ancestorNodes.addClass('highlighted')
      ancestorEdges.addClass('highlighted')

      const goNode = nodeMap.get(node.id())
      if (goNode) {
        const renderedPos = node.renderedPosition()
        setTooltip({ x: renderedPos.x, y: renderedPos.y, node: goNode })
      }
    })

    cy.on('mouseout', 'node', () => {
      cy.nodes().removeClass('highlighted')
      cy.edges().removeClass('highlighted')
      setTooltip(null)
    })

    // Click: notify parent for peptide pane selection
    cy.on('tap', 'node', (evt) => {
      onNodeClickRef.current?.(evt.target.id())
    })

    // Click background: clear selection
    cy.on('tap', (evt) => {
      if (evt.target === cy) {
        onNodeClickRef.current?.(null)
      }
    })

    // Dismiss tooltip on zoom/pan
    cy.on('zoom pan', () => {
      setTooltip(null)
    })

    cyRef.current = cy

    return () => {
      cy.destroy()
      cyRef.current = null
    }
  }, [nodes, buildElements, baseColor, isDark])

  // Update colors when metric or baseColor changes (without re-layout)
  useEffect(() => {
    const cy = cyRef.current
    if (!cy || nodes.length === 0) return

    const normalized = normalizeValues(nodes, metric)
    const stops = generateColorStops(baseColor, isDark)
    cy.batch(() => {
      cy.nodes().forEach((ele) => {
        const t = normalized.get(ele.id()) ?? 0
        const bgColor = interpolateColor(t, stops)
        ele.style('background-color', bgColor)
        ele.style('color', textColorForBg(bgColor, isDark))
      })
    })
  }, [metric, nodes, baseColor, isDark])

  // Export PNG handler — exposed via ref or callback
  useEffect(() => {
    // Attach export function to the container element for parent access
    const el = containerRef.current
    const exportBg = (isDark ? GO_DAG.dark : GO_DAG.light).exportBg
    if (el) {
      (el as any).__exportPng = () => {
        const cy = cyRef.current
        if (!cy) return
        const png = cy.png({ full: true, scale: 2, bg: exportBg })
        const link = document.createElement('a')
        link.href = png
        link.download = 'go_dag.png'
        link.click()
      }
      ;(el as any).__exportSvg = () => {
        const cy = cyRef.current
        if (!cy) return
        const svgContent = (cy as any).svg({ full: true, scale: 1, bg: exportBg })
        const blob = new Blob([svgContent], { type: 'image/svg+xml;charset=utf-8' })
        const url = URL.createObjectURL(blob)
        const link = document.createElement('a')
        link.href = url
        link.download = 'go_dag.svg'
        link.click()
        URL.revokeObjectURL(url)
      }
    }
  })

  if (nodes.length === 0) {
    return (
      <div className="flex items-center justify-center h-full text-gray-500 dark:text-gray-400">
        No GO terms in this namespace.
      </div>
    )
  }

  return (
    <div className="relative w-full h-full overflow-visible">
      <div ref={containerRef} data-cy-container className="absolute inset-0" />

      {/* Tooltip */}
      {tooltip && (
        <div
          ref={tooltipRef}
          className="fixed z-50 bg-white dark:bg-gray-800 border border-gray-300 dark:border-gray-600 rounded-lg shadow-lg dark:shadow-indigo-500/10 p-3 text-sm pointer-events-none"
          style={{
            left: tooltipPos.left,
            top: tooltipPos.top,
            maxWidth: 320,
          }}
        >
          <p className="font-mono text-xs text-gray-500 dark:text-gray-400 mb-1">{tooltip.node.id}</p>
          <p className="font-semibold text-gray-900 dark:text-gray-100 mb-2">{tooltip.node.name}</p>
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
            {filterLabel && tooltip.node.qvalueGoForTaxon != null && (
              <>
                <span className="text-gray-600 dark:text-gray-400">Q-value (GO for Taxon):</span>
                <span className="font-medium text-right text-gray-900 dark:text-gray-100">{tooltip.node.qvalueGoForTaxon < 0.001 ? tooltip.node.qvalueGoForTaxon.toExponential(3) : tooltip.node.qvalueGoForTaxon.toFixed(4)}</span>
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
