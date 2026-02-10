import { useCallback, useEffect, useMemo, useRef, useState } from 'react'
import { Link, useParams } from 'react-router-dom'
import { ArrowLeft, Loader2, AlertCircle } from 'lucide-react'
import TaxonomyChart from '../components/TaxonomyChart'
import TaxonomyControls from '../components/TaxonomyControls'
import { parseTaxonomyCsv, filterCanonicalRanks, filterByMaxRank, validateCanonicalHierarchy, ensureStrictRankLayers } from '../utils/taxonomyParser'
import type { TaxonNode, CanonicalRank } from '../utils/taxonomyParser'

export type { TaxonNode } from '../utils/taxonomyParser'
export type ChartType = 'sunburst' | 'treemap' | 'icicle' | 'sankey'

export default function TaxonomyPage() {
  const { jobId, listId } = useParams<{ jobId: string; listId: string }>()
  const [allNodes, setAllNodes] = useState<TaxonNode[]>([])
  const [canonicalNodes, setCanonicalNodes] = useState<TaxonNode[]>([])
  const [loading, setLoading] = useState(true)
  const [error, setError] = useState<string | null>(null)

  const [chartType, setChartType] = useState<ChartType>('sunburst')
  const [maxRank, setMaxRank] = useState<CanonicalRank>('species')
  const [minRatioTotal, setMinRatioTotal] = useState(0.001)
  const chartContainerRef = useRef<HTMLDivElement>(null)

  const filteredNodes = useMemo(
    () => {
      // Filter by min ratio, but always keep the root node
      const byRatio = canonicalNodes.filter(n => n.rank === 'root' || n.ratioTotal >= minRatioTotal)
      const byRank = filterByMaxRank(byRatio, maxRank)
      // Sankey uses explicit x positions for rank alignment — no placeholders needed
      return chartType === 'sankey' ? byRank : ensureStrictRankLayers(byRank)
    },
    [canonicalNodes, maxRank, minRatioTotal, chartType],
  )

  const handleExportPng = useCallback(() => {
    const el = chartContainerRef.current?.querySelector('.js-plotly-plot') as any
    if (el) {
      import('plotly.js-dist-min').then((Plotly: any) => {
        Plotly.default.downloadImage(el, {
          format: 'png',
          width: 1600,
          height: 1200,
          filename: 'taxonomy',
          scale: 2,
        })
      })
    }
  }, [])

  const handleExportSvg = useCallback(() => {
    const el = chartContainerRef.current?.querySelector('.js-plotly-plot') as any
    if (el) {
      import('plotly.js-dist-min').then((Plotly: any) => {
        Plotly.default.downloadImage(el, {
          format: 'svg',
          width: 1600,
          height: 1200,
          filename: 'taxonomy',
        })
      })
    }
  }, [])

  useEffect(() => {
    if (!jobId || !listId) return

    async function fetchData() {
      try {
        const response = await fetch(`/api/jobs/${jobId}/results/${listId}/taxonomy_nodes.csv`)
        if (!response.ok) throw new Error(`Failed to fetch taxonomy data: ${response.status}`)
        const text = await response.text()
        const nodes = parseTaxonomyCsv(text)
        setAllNodes(nodes)
        const canonical = filterCanonicalRanks(nodes)
        const validationErrors = validateCanonicalHierarchy(canonical)
        if (validationErrors.length > 0) {
          throw new Error(
            `Taxonomy hierarchy validation failed:\n${validationErrors.join('\n')}`
          )
        }
        setCanonicalNodes(canonical)
        setError(null)
      } catch (err) {
        setError(err instanceof Error ? err.message : 'Unknown error')
      } finally {
        setLoading(false)
      }
    }

    fetchData()
  }, [jobId, listId])

  if (loading) {
    return (
      <div className="flex items-center justify-center py-12">
        <Loader2 className="w-8 h-8 text-indigo-600 animate-spin" />
      </div>
    )
  }

  if (error) {
    return (
      <div className="text-center py-12">
        <AlertCircle className="w-12 h-12 text-red-500 mx-auto mb-4" />
        <h2 className="text-xl font-semibold text-gray-900 mb-2">Error</h2>
        <p className="text-gray-600">{error}</p>
        <Link to={`/job/${jobId}`} className="mt-4 inline-flex items-center gap-1 text-indigo-600 hover:text-indigo-800">
          <ArrowLeft className="w-4 h-4" /> Back to Job
        </Link>
      </div>
    )
  }

  return (
    <div className="flex flex-col h-[calc(100vh-8rem)]">
      {/* Header */}
      <div className="flex items-center justify-between mb-4">
        <div className="flex items-center gap-4">
          <Link
            to={`/job/${jobId}`}
            className="inline-flex items-center gap-1 text-sm text-gray-600 hover:text-gray-900 transition-colors"
          >
            <ArrowLeft className="w-4 h-4" />
            Back to Job
          </Link>
          <h1 className="text-xl font-bold text-gray-900">
            Taxonomy Visualization
          </h1>
          <span className="text-sm text-gray-500">
            {listId} · {canonicalNodes.length} terms ({allNodes.length} total)
          </span>
        </div>
      </div>

      {/* Controls */}
      <TaxonomyControls
        chartType={chartType}
        onChartTypeChange={setChartType}
        maxRank={maxRank}
        onMaxRankChange={setMaxRank}
        minRatioTotal={minRatioTotal}
        onMinRatioTotalChange={setMinRatioTotal}
        filteredNodeCount={filteredNodes.length}
        totalNodeCount={canonicalNodes.length}
        onExportPng={handleExportPng}
        onExportSvg={handleExportSvg}
      />

      {/* Chart */}
      <div className="flex-1 min-h-0 border border-gray-200 rounded-lg overflow-hidden bg-white">
        <div ref={chartContainerRef} className="w-full h-full">
          <TaxonomyChart
            nodes={filteredNodes}
            chartType={chartType}
          />
        </div>
      </div>
    </div>
  )
}
