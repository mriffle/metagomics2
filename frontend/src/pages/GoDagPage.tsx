import { useCallback, useEffect, useRef, useState } from 'react'
import { Link, useParams } from 'react-router-dom'
import { ArrowLeft, Loader2, AlertCircle } from 'lucide-react'
import GoDagViewer from '../components/GoDagViewer'
import GoDagControls from '../components/GoDagControls'
import { parseGoTermsCsv } from '../utils/goParser'
import type { GoTermNode } from '../utils/goParser'

export type { GoTermNode } from '../utils/goParser'
export type MetricKey = 'quantity' | 'ratioTotal' | 'ratioAnnotated' | 'nPeptides'

const NAMESPACE_LABELS: Record<string, string> = {
  biological_process: 'Biological Process',
  cellular_component: 'Cellular Component',
  molecular_function: 'Molecular Function',
}

export default function GoDagPage() {
  const { jobId, listId } = useParams<{ jobId: string; listId: string }>()
  const [allNodes, setAllNodes] = useState<GoTermNode[]>([])
  const [loading, setLoading] = useState(true)
  const [error, setError] = useState<string | null>(null)
  const [danglingParents, setDanglingParents] = useState<string[]>([])
  const [listFilename, setListFilename] = useState<string>('')

  const [selectedNamespace, setSelectedNamespace] = useState('biological_process')
  const [selectedMetric, setSelectedMetric] = useState<MetricKey>('quantity')
  const [minRatioTotal, setMinRatioTotal] = useState(0.1)
  const viewerContainerRef = useRef<HTMLDivElement>(null)

  const handleExportPng = useCallback(() => {
    const el = viewerContainerRef.current?.querySelector('[data-cy-container]') as any
    if (el?.__exportPng) el.__exportPng()
  }, [])

  const handleExportSvg = useCallback(() => {
    const el = viewerContainerRef.current?.querySelector('[data-cy-container]') as any
    if (el?.__exportSvg) el.__exportSvg()
  }, [])

  useEffect(() => {
    if (!jobId || !listId) return

    // Fetch the job to get the filename for this list
    fetch(`/api/jobs/${jobId}`)
      .then(res => res.ok ? res.json() : null)
      .then(data => {
        if (data?.peptide_lists) {
          const match = data.peptide_lists.find((pl: any) => pl.list_id === listId)
          if (match) setListFilename(match.filename)
        }
      })
      .catch(() => {})

    async function fetchData() {
      try {
        const response = await fetch(`/api/jobs/${jobId}/results/${listId}/go_terms.csv`)
        if (!response.ok) throw new Error(`Failed to fetch GO terms: ${response.status}`)
        const text = await response.text()
        const nodes = parseGoTermsCsv(text)
        setAllNodes(nodes)

        // Check for dangling parents
        const nodeIds = new Set(nodes.map(n => n.id))
        const dangling = new Set<string>()
        for (const node of nodes) {
          for (const pid of node.parentIds) {
            if (!nodeIds.has(pid)) dangling.add(pid)
          }
        }
        setDanglingParents(Array.from(dangling).sort())
        setError(null)
      } catch (err) {
        setError(err instanceof Error ? err.message : 'Unknown error')
      } finally {
        setLoading(false)
      }
    }

    fetchData()
  }, [jobId, listId])

  // Get available namespaces from data
  const namespaces = Array.from(new Set(allNodes.map(n => n.namespace))).filter(Boolean).sort()

  // Auto-select first available namespace if current selection has no data
  useEffect(() => {
    if (namespaces.length > 0 && !namespaces.includes(selectedNamespace)) {
      setSelectedNamespace(namespaces[0])
    }
  }, [namespaces])

  // Filter nodes by selected namespace and abundance cutoff
  const filteredNodes = allNodes.filter(
    n => n.namespace === selectedNamespace && n.ratioTotal >= minRatioTotal
  )

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
            Gene Ontology DAG
          </h1>
          <span className="text-sm text-gray-500">
            {listFilename || listId} · {allNodes.length} terms
          </span>
        </div>
      </div>

      {/* Dangling parent warning */}
      {danglingParents.length > 0 && (
        <div className="mb-3 p-3 bg-amber-50 border border-amber-200 rounded-lg text-sm text-amber-800">
          <span className="font-medium">Data integrity warning:</span>{' '}
          {danglingParents.length} parent GO ID(s) referenced but not present in the data:{' '}
          <span className="font-mono text-xs">{danglingParents.slice(0, 10).join(', ')}</span>
          {danglingParents.length > 10 && <span> and {danglingParents.length - 10} more</span>}
        </div>
      )}

      {/* Controls */}
      <GoDagControls
        namespaces={namespaces}
        namespaceLabels={NAMESPACE_LABELS}
        selectedNamespace={selectedNamespace}
        onNamespaceChange={setSelectedNamespace}
        selectedMetric={selectedMetric}
        onMetricChange={setSelectedMetric}
        minRatioTotal={minRatioTotal}
        onMinRatioTotalChange={setMinRatioTotal}
        filteredNodeCount={filteredNodes.length}
        totalNodeCount={allNodes.filter(n => n.namespace === selectedNamespace).length}
        onExportPng={handleExportPng}
        onExportSvg={handleExportSvg}
      />

      {/* Graph */}
      <div className="flex-1 min-h-0 border border-gray-200 rounded-lg overflow-hidden bg-white">
        <div ref={viewerContainerRef} className="w-full h-full">
          <GoDagViewer
            nodes={filteredNodes}
            metric={selectedMetric}
            key={selectedNamespace}
          />
        </div>
      </div>
    </div>
  )
}
