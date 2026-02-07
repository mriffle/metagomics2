import { useEffect, useRef, useState } from 'react'
import { useParams } from 'react-router-dom'
import {
  Clock,
  CheckCircle,
  XCircle,
  Loader2,
  Download,
  FileText,
  AlertCircle,
} from 'lucide-react'

interface PeptideList {
  list_id: string
  filename: string
  status: string
  n_peptides: number | null
  n_matched: number | null
  n_unmatched: number | null
}

interface Job {
  job_id: string
  created_at: string
  status: string
  progress_done: number
  progress_total: number
  current_step: string | null
  error_message: string | null
  peptide_lists: PeptideList[]
}

const statusIcons: Record<string, React.ReactNode> = {
  uploaded: <Clock className="w-6 h-6 text-gray-400" />,
  queued: <Clock className="w-6 h-6 text-yellow-500" />,
  running: <Loader2 className="w-6 h-6 text-blue-500 animate-spin" />,
  completed: <CheckCircle className="w-6 h-6 text-green-500" />,
  failed: <XCircle className="w-6 h-6 text-red-500" />,
}

export default function JobPage() {
  const { jobId } = useParams<{ jobId: string }>()
  const [job, setJob] = useState<Job | null>(null)
  const [loading, setLoading] = useState(true)
  const [error, setError] = useState<string | null>(null)

  const intervalRef = useRef<ReturnType<typeof setInterval> | null>(null)

  useEffect(() => {
    if (!jobId) return

    fetchJob()
    intervalRef.current = setInterval(fetchJob, 3000)
    return () => {
      if (intervalRef.current) clearInterval(intervalRef.current)
    }
  }, [jobId])

  // Stop polling once job reaches a terminal state
  useEffect(() => {
    if (job && (job.status === 'completed' || job.status === 'failed')) {
      if (intervalRef.current) {
        clearInterval(intervalRef.current)
        intervalRef.current = null
      }
    }
  }, [job?.status])

  async function fetchJob() {
    try {
      const response = await fetch(`/api/jobs/${jobId}`)
      if (!response.ok) {
        if (response.status === 404) {
          throw new Error('Job not found')
        }
        throw new Error('Failed to fetch job')
      }
      const data = await response.json()
      setJob(data)
      setError(null)
    } catch (err) {
      setError(err instanceof Error ? err.message : 'Unknown error')
    } finally {
      setLoading(false)
    }
  }

  if (loading) {
    return (
      <div className="flex items-center justify-center py-12">
        <Loader2 className="w-8 h-8 text-indigo-600 animate-spin" />
      </div>
    )
  }

  if (error || !job) {
    return (
      <div className="text-center py-12">
        <AlertCircle className="w-12 h-12 text-red-500 mx-auto mb-4" />
        <h2 className="text-xl font-semibold text-gray-900 mb-2">Error</h2>
        <p className="text-gray-600">{error || 'Job not found'}</p>
      </div>
    )
  }

  const progressPercent = job.progress_total > 0
    ? Math.round((job.progress_done / job.progress_total) * 100)
    : 0

  return (
    <div className="space-y-8">
      {/* Header */}
      <div className="bg-white rounded-lg border border-gray-200 p-6">
        <div className="flex items-start justify-between">
          <div>
            <h1 className="text-2xl font-bold text-gray-900 flex items-center gap-3">
              {statusIcons[job.status]}
              Job Status
            </h1>
            <p className="mt-1 text-sm text-gray-500 font-mono">{job.job_id}</p>
          </div>
          <div className="text-right">
            <p className="text-sm text-gray-500">Created</p>
            <p className="text-gray-900">{new Date(job.created_at).toLocaleString()}</p>
          </div>
        </div>

        {/* Progress Bar */}
        {(job.status === 'running' || job.status === 'queued') && (
          <div className="mt-6">
            <div className="flex justify-between text-sm text-gray-600 mb-2">
              <span>{job.current_step || 'Processing...'}</span>
              <span>{job.progress_done} / {job.progress_total} peptide lists</span>
            </div>
            <div className="w-full bg-gray-200 rounded-full h-3">
              <div
                className="bg-indigo-600 h-3 rounded-full transition-all duration-500"
                style={{ width: `${progressPercent}%` }}
              />
            </div>
          </div>
        )}

        {/* Error Message */}
        {job.status === 'failed' && job.error_message && (
          <div className="mt-6 p-4 bg-red-50 border border-red-200 rounded-lg text-red-700">
            <p className="font-medium">Error:</p>
            <p className="mt-1">{job.error_message}</p>
          </div>
        )}
      </div>

      {/* Peptide Lists */}
      <div className="bg-white rounded-lg border border-gray-200 p-6">
        <h2 className="text-lg font-semibold text-gray-900 mb-4 flex items-center gap-2">
          <FileText className="w-5 h-5" />
          Peptide Lists
        </h2>

        <div className="space-y-4">
          {job.peptide_lists.map((list) => (
            <div
              key={list.list_id}
              className="border border-gray-200 rounded-lg p-4"
            >
              <div className="flex items-center justify-between">
                <div>
                  <p className="font-medium text-gray-900">{list.filename}</p>
                  <p className="text-sm text-gray-500">{list.list_id}</p>
                </div>
                <div className="flex items-center gap-4">
                  {list.n_peptides !== null && (
                    <div className="text-sm text-gray-600">
                      <span className="font-medium">{list.n_peptides}</span> peptides
                      {list.n_matched !== null && (
                        <span className="ml-2">
                          ({list.n_matched} matched)
                        </span>
                      )}
                    </div>
                  )}
                  <span className={`px-2 py-1 rounded text-xs font-medium ${
                    list.status === 'done' ? 'bg-green-100 text-green-800' :
                    list.status === 'running' ? 'bg-blue-100 text-blue-800' :
                    list.status === 'failed' ? 'bg-red-100 text-red-800' :
                    'bg-gray-100 text-gray-800'
                  }`}>
                    {list.status}
                  </span>
                </div>
              </div>
            </div>
          ))}
        </div>
      </div>

      {/* Downloads */}
      {job.status === 'completed' && (
        <div className="bg-white rounded-lg border border-gray-200 p-6">
          <h2 className="text-lg font-semibold text-gray-900 mb-4 flex items-center gap-2">
            <Download className="w-5 h-5" />
            Downloads
          </h2>

          <div className="mb-4 text-sm text-gray-600 space-y-1">
            <p><span className="font-medium">taxonomy_nodes.csv</span> — Taxonomy assignments with aggregated peptide quantities and ratios at each node.</p>
            <p><span className="font-medium">go_terms.csv</span> — Gene Ontology term assignments with aggregated peptide quantities and ratios.</p>
            <p><span className="font-medium">coverage.csv</span> — Summary of how many peptides (and what fraction of total quantity) received annotation.</p>
            <p><span className="font-medium">run_manifest.json</span> — Provenance record of parameters, software version, and reference data used for this run.</p>
          </div>

          <div className="space-y-4">
            {/* Download All */}
            <a
              href={`/api/jobs/${job.job_id}/results/all_results.zip`}
              className="block p-4 border border-gray-200 rounded-lg hover:bg-gray-50 transition-colors"
            >
              <div className="flex items-center gap-3">
                <Download className="w-8 h-8 text-indigo-600" />
                <div>
                  <p className="font-medium text-gray-900">Download All Results</p>
                  <p className="text-sm text-gray-500">ZIP archive with all output files</p>
                </div>
              </div>
            </a>

            {/* Individual Files */}
            <div className="grid grid-cols-2 gap-4">
              {job.peptide_lists.map((list) => (
                <div key={list.list_id} className="space-y-2">
                  <p className="text-sm font-medium text-gray-700">{list.filename}</p>
                  <div className="flex flex-wrap gap-2">
                    {['taxonomy_nodes.csv', 'go_terms.csv', 'coverage.csv', 'run_manifest.json'].map((file) => (
                      <a
                        key={file}
                        href={`/api/jobs/${job.job_id}/results/${list.list_id}/${file}`}
                        className="inline-flex items-center gap-1 px-2 py-1 text-xs bg-gray-100 hover:bg-gray-200 rounded transition-colors"
                      >
                        <FileText className="w-3 h-3" />
                        {file}
                      </a>
                    ))}
                  </div>
                </div>
              ))}
            </div>
          </div>
        </div>
      )}
    </div>
  )
}
