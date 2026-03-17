import { useEffect, useRef, useState } from 'react'
import { Link, useNavigate, useParams } from 'react-router-dom'
import {
  Clock,
  CheckCircle,
  XCircle,
  Loader2,
  Download,
  FileText,
  AlertCircle,
  GitBranch,
  TreePine,
  RefreshCw,
} from 'lucide-react'
import { STATUS_BADGE_CLASSES, STATUS_BADGE_DEFAULT } from '../utils/colors'

interface PeptideList {
  list_id: string
  filename: string
  status: string
  n_peptides: number | null
  n_matched: number | null
  n_unmatched: number | null
}

interface JobParams {
  search_tool: string
  db_choice: string
  db_name: string
  max_evalue: number | null
  min_pident: number | null
  min_qcov: number | null
  min_alnlen: number | null
  top_k: number | null
  notification_email: string
  fasta_filename: string
}

interface Job {
  job_id: string
  created_at: string
  status: string
  params: JobParams
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
  const navigate = useNavigate()
  const [job, setJob] = useState<Job | null>(null)
  const [loading, setLoading] = useState(true)
  const [error, setError] = useState<string | null>(null)
  const [regenerating, setRegenerating] = useState(false)

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

  async function handleRegenerateId() {
    if (!jobId || regenerating) return
    if (!window.confirm('This will change the URL for this job. The old URL will stop working. Continue?')) return
    setRegenerating(true)
    try {
      const response = await fetch(`/api/jobs/${jobId}/regenerate-id`, { method: 'POST' })
      if (!response.ok) {
        const data = await response.json()
        throw new Error(data.detail || 'Failed to regenerate ID')
      }
      const data = await response.json()
      navigate(`/job/${data.new_job_id}`, { replace: true })
    } catch (err) {
      setError(err instanceof Error ? err.message : 'Unknown error')
    } finally {
      setRegenerating(false)
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
        <h2 className="text-xl font-semibold text-gray-900 dark:text-gray-100 mb-2">Error</h2>
        <p className="text-gray-600 dark:text-gray-400">{error || 'Job not found'}</p>
      </div>
    )
  }

  const progressPercent = job.progress_total > 0
    ? Math.round((job.progress_done / job.progress_total) * 100)
    : 0

  return (
    <div className="space-y-8">
      {/* Header */}
      <div className="bg-white dark:bg-gray-900 rounded-lg border border-gray-200 dark:border-gray-700 p-6">
        <div className="flex items-start justify-between">
          <div>
            <h1 className="text-2xl font-bold text-gray-900 dark:text-gray-100 flex items-center gap-3">
              {statusIcons[job.status]}
              Job Status
            </h1>
            <div className="mt-1 flex items-center gap-2">
              <p className="text-sm text-gray-500 dark:text-gray-400 font-mono">{job.job_id}</p>
              <span className="relative group inline-flex">
                <button
                  onClick={handleRegenerateId}
                  disabled={regenerating}
                  className="inline-flex items-center gap-1 px-2 py-0.5 text-xs bg-gray-100 dark:bg-gray-800 hover:bg-gray-200 dark:hover:bg-gray-700 text-gray-600 dark:text-gray-300 rounded transition-colors disabled:opacity-50"
                >
                  <RefreshCw className={`w-3 h-3 ${regenerating ? 'animate-spin' : ''}`} />
                  Change Hash
                </button>
                <span className="absolute bottom-full left-1/2 -translate-x-1/2 mb-2 px-3 py-2 text-xs text-white bg-gray-800 dark:bg-gray-700 rounded-lg whitespace-normal w-56 text-center opacity-0 invisible group-hover:opacity-100 group-hover:visible transition-all z-10 pointer-events-none">
                  Generate a new URL for this job. The old URL will stop working. Use this to revoke access if you previously shared the link.
                  <span className="absolute top-full left-1/2 -translate-x-1/2 border-4 border-transparent border-t-gray-800 dark:border-t-gray-700" />
                </span>
              </span>
            </div>
          </div>
          <div className="text-right">
            <p className="text-sm text-gray-500 dark:text-gray-400">Created</p>
            <p className="text-gray-900 dark:text-gray-100">{new Date(job.created_at).toLocaleString()}</p>
          </div>
        </div>

        {/* Parameters */}
        <div className="mt-4 flex flex-wrap gap-x-6 gap-y-1 text-sm text-gray-600 dark:text-gray-400">
          {job.params.fasta_filename && (
            <span><span className="font-medium text-gray-700 dark:text-gray-300">FASTA:</span> {job.params.fasta_filename}</span>
          )}
          {job.params.db_choice && (
            <span><span className="font-medium text-gray-700 dark:text-gray-300">Database:</span> {job.params.db_name ? `${job.params.db_name} (${job.params.db_choice})` : job.params.db_choice}</span>
          )}
          {job.params.max_evalue != null && (
            <span><span className="font-medium text-gray-700 dark:text-gray-300">Max E-value:</span> {job.params.max_evalue}</span>
          )}
          {job.params.min_pident != null && (
            <span><span className="font-medium text-gray-700 dark:text-gray-300">Min % Identity:</span> {job.params.min_pident}</span>
          )}
          {job.params.top_k != null && (
            <span><span className="font-medium text-gray-700 dark:text-gray-300">Top K:</span> {job.params.top_k}</span>
          )}
          {job.params.min_qcov != null && (
            <span><span className="font-medium text-gray-700 dark:text-gray-300">Min Query Cov:</span> {job.params.min_qcov}</span>
          )}
          {job.params.min_alnlen != null && (
            <span><span className="font-medium text-gray-700 dark:text-gray-300">Min Aln Length:</span> {job.params.min_alnlen}</span>
          )}
        </div>

        {/* Progress Bar */}
        {(job.status === 'running' || job.status === 'queued') && (
          <div className="mt-6">
            <div className="flex justify-between text-sm text-gray-600 dark:text-gray-400 mb-2">
              <span>{job.current_step || 'Processing...'}</span>
              <span>{progressPercent}%</span>
            </div>
            <div className="w-full bg-gray-200 dark:bg-gray-700 rounded-full h-3">
              <div
                className="bg-indigo-600 h-3 rounded-full transition-all duration-500"
                style={{ width: `${progressPercent}%` }}
              />
            </div>
          </div>
        )}

        {/* Error Message */}
        {job.status === 'failed' && job.error_message && (
          <div className="mt-6 p-4 bg-red-50 dark:bg-red-900/30 border border-red-200 dark:border-red-800 rounded-lg text-red-700 dark:text-red-400">
            <p className="font-medium">Error:</p>
            <p className="mt-1">{job.error_message}</p>
          </div>
        )}
      </div>

      {/* Peptide Lists */}
      <div className="bg-white dark:bg-gray-900 rounded-lg border border-gray-200 dark:border-gray-700 p-6">
        <div className="flex items-center justify-between mb-4">
          <h2 className="text-lg font-semibold text-gray-900 dark:text-gray-100 flex items-center gap-2">
            <FileText className="w-5 h-5" />
            Peptide Lists
          </h2>
          {job.status === 'completed' && (
            <a
              href={`/api/jobs/${job.job_id}/results/all_results.zip`}
              className="inline-flex items-center gap-2 px-3 py-1.5 text-sm bg-indigo-600 hover:bg-indigo-700 text-white rounded-md transition-colors"
            >
              <Download className="w-4 h-4" />
              Download All Results (ZIP)
            </a>
          )}
        </div>

        <div className="space-y-4">
          {job.peptide_lists.map((list) => (
            <div
              key={list.list_id}
              className="border border-gray-200 dark:border-gray-700 rounded-lg p-4"
            >
              <div className="flex items-center justify-between">
                <div>
                  <p className="font-medium text-gray-900 dark:text-gray-100">{list.filename}</p>
                </div>
                <div className="flex items-center gap-4">
                  {list.n_peptides !== null && (
                    <div className="text-sm text-gray-600 dark:text-gray-400">
                      <span className="font-medium">{list.n_peptides}</span> peptides
                      {list.n_matched !== null && (
                        <span className="ml-2">
                          ({list.n_matched} matched)
                        </span>
                      )}
                    </div>
                  )}
                  <span className={`px-2 py-1 rounded text-xs font-medium ${
                    STATUS_BADGE_CLASSES[list.status] ?? STATUS_BADGE_DEFAULT
                  }`}>
                    {list.status}
                  </span>
                </div>
              </div>

              {/* Per-list downloads and visualizations */}
              {job.status === 'completed' && list.status === 'done' && (
                <div className="mt-3 pt-3 border-t border-gray-100 dark:border-gray-700 space-y-2">
                  <div className="flex items-center gap-2 flex-wrap">
                    <span className="text-xs font-medium text-gray-500 dark:text-gray-400">Download results:</span>
                    {[
                      { name: 'taxonomy_nodes.csv', tip: 'Taxonomy assignments with aggregated peptide quantities and ratios at each node' },
                      { name: 'go_terms.csv', tip: 'Gene Ontology term assignments with aggregated peptide quantities and ratios' },
                      { name: 'go_taxonomy_combo.csv', tip: 'Cross-tabulation of taxonomy nodes and GO terms with quantities and fractions for each combination' },
                      { name: 'coverage.csv', tip: 'Summary of how many peptides (and what fraction of total quantity) received any annotations (taxonomy or GO)' },
                      { name: 'run_manifest.json', tip: 'Provenance record of parameters, software versions, and reference data used for this run' },
                      { name: 'peptide_mapping.parquet', tip: 'Full peptide-to-protein mapping table with taxonomy and GO annotations. Query with DuckDB or any Parquet-compatible tool.' },
                    ].map((file) => (
                      <span key={file.name} className="relative group inline-flex">
                        <a
                          href={`/api/jobs/${job.job_id}/results/${list.list_id}/${file.name}`}
                          className="inline-flex items-center gap-1 px-2 py-1 text-xs bg-gray-100 dark:bg-gray-800 hover:bg-gray-200 dark:hover:bg-gray-700 text-gray-700 dark:text-gray-300 rounded transition-colors"
                        >
                          <FileText className="w-3 h-3" />
                          {file.name}
                        </a>
                        <span className="absolute bottom-full left-1/2 -translate-x-1/2 mb-2 px-3 py-2 text-xs text-white bg-gray-800 dark:bg-gray-700 rounded-lg whitespace-normal w-56 text-center opacity-0 invisible group-hover:opacity-100 group-hover:visible transition-all z-10 pointer-events-none">
                          {file.tip}
                          <span className="absolute top-full left-1/2 -translate-x-1/2 border-4 border-transparent border-t-gray-800 dark:border-t-gray-700" />
                        </span>
                      </span>
                    ))}
                  </div>
                  <div className="flex items-center gap-2 flex-wrap">
                    <span className="text-xs font-medium text-gray-500 dark:text-gray-400">View results:</span>
                    <span className="relative group inline-flex">
                      <Link
                        to={`/job/${job.job_id}/go/${list.list_id}`}
                        className="inline-flex items-center gap-1 px-2 py-1 text-xs bg-indigo-100 dark:bg-indigo-900/50 hover:bg-indigo-200 dark:hover:bg-indigo-800/50 text-indigo-700 dark:text-indigo-300 rounded transition-colors"
                      >
                        <GitBranch className="w-3 h-3" />
                        Gene Ontology
                      </Link>
                      <span className="absolute bottom-full left-1/2 -translate-x-1/2 mb-2 px-3 py-2 text-xs text-white bg-gray-800 dark:bg-gray-700 rounded-lg whitespace-normal w-56 text-center opacity-0 invisible group-hover:opacity-100 group-hover:visible transition-all z-10 pointer-events-none">
                        Visualize Gene Ontology terms found for this peptide list
                        <span className="absolute top-full left-1/2 -translate-x-1/2 border-4 border-transparent border-t-gray-800 dark:border-t-gray-700" />
                      </span>
                    </span>
                    <span className="relative group inline-flex">
                      <Link
                        to={`/job/${job.job_id}/taxonomy/${list.list_id}`}
                        className="inline-flex items-center gap-1 px-2 py-1 text-xs bg-emerald-100 dark:bg-emerald-900/50 hover:bg-emerald-200 dark:hover:bg-emerald-800/50 text-emerald-700 dark:text-emerald-300 rounded transition-colors"
                      >
                        <TreePine className="w-3 h-3" />
                        Taxonomy
                      </Link>
                      <span className="absolute bottom-full left-1/2 -translate-x-1/2 mb-2 px-3 py-2 text-xs text-white bg-gray-800 dark:bg-gray-700 rounded-lg whitespace-normal w-56 text-center opacity-0 invisible group-hover:opacity-100 group-hover:visible transition-all z-10 pointer-events-none">
                        Visualize taxonomic assignments found for this peptide list
                        <span className="absolute top-full left-1/2 -translate-x-1/2 border-4 border-transparent border-t-gray-800 dark:border-t-gray-700" />
                      </span>
                    </span>
                  </div>
                </div>
              )}
            </div>
          ))}
        </div>
      </div>
    </div>
  )
}
