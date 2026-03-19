import { useEffect, useState } from 'react'
import { useNavigate } from 'react-router-dom'
import { Upload, FileText, Settings, Loader2, AlertCircle, HelpCircle, FlaskConical } from 'lucide-react'

function Tooltip({ text }: { text: string }) {
  return (
    <span className="relative group inline-flex ml-1">
      <HelpCircle className="w-4 h-4 text-gray-400 dark:text-gray-500 cursor-help" />
      <span className="absolute bottom-full left-1/2 -translate-x-1/2 mb-2 px-3 py-2 text-xs text-white bg-gray-800 dark:bg-gray-700 rounded-lg whitespace-normal w-56 text-center opacity-0 invisible group-hover:opacity-100 group-hover:visible transition-all z-10 pointer-events-none">
        {text}
        <span className="absolute top-full left-1/2 -translate-x-1/2 border-4 border-transparent border-t-gray-800 dark:border-t-gray-700" />
      </span>
    </span>
  )
}

export default function NewJobPage() {
  const navigate = useNavigate()
  const [fastaFile, setFastaFile] = useState<File | null>(null)
  const [peptideFiles, setPeptideFiles] = useState<File[]>([])
  const [submitting, setSubmitting] = useState(false)
  const [error, setError] = useState<string | null>(null)

  // Config from server
  const [diamondVersion, setDiamondVersion] = useState('')
  const [databases, setDatabases] = useState<{ name: string; description: string; path: string }[]>([])

  useEffect(() => {
    fetch('/api/config')
      .then(res => res.json())
      .then(data => {
        setDiamondVersion(data.diamond_version || '')
        const dbs = data.databases || []
        setDatabases(dbs)
        if (dbs.length > 0) {
          setDbChoice(dbs[0].path)
        }
      })
      .catch(() => {
        setDiamondVersion('')
        setDatabases([])
      })
  }, [])

  // Parameters
  const [searchTool, setSearchTool] = useState('diamond')
  const [dbChoice, setDbChoice] = useState('')
  const [maxEvalue, setMaxEvalue] = useState('1e-10')
  const [minPident, setMinPident] = useState('80')
  const [topK, setTopK] = useState('1')
  const [notificationEmail, setNotificationEmail] = useState('')
  const [computeEnrichment, setComputeEnrichment] = useState(false)
  const [enrichmentIterations, setEnrichmentIterations] = useState('1000')

  async function handleSubmit(e: React.FormEvent) {
    e.preventDefault()

    if (!fastaFile) {
      setError('Please select a FASTA file')
      return
    }
    if (peptideFiles.length === 0) {
      setError('Please select at least one peptide file')
      return
    }
    setSubmitting(true)
    setError(null)

    try {
      const formData = new FormData()
      formData.append('fasta', fastaFile)
      peptideFiles.forEach((file) => {
        formData.append('peptides', file)
      })

      const params: Record<string, unknown> = {
        search_tool: searchTool,
        db_choice: dbChoice,
      }
      if (maxEvalue) params.max_evalue = parseFloat(maxEvalue)
      if (minPident) params.min_pident = parseFloat(minPident)
      if (topK) params.top_k = parseInt(topK)
      const emailTrimmed = notificationEmail.trim()
      if (emailTrimmed) params.notification_email = emailTrimmed
      if (computeEnrichment) {
        params.compute_enrichment_pvalues = true
        const iters = parseInt(enrichmentIterations)
        if (!isNaN(iters) && iters >= 100 && iters <= 10000) {
          params.enrichment_iterations = iters
        }
      }

      formData.append('params', JSON.stringify(params))

      const response = await fetch('/api/jobs', {
        method: 'POST',
        body: formData,
      })

      if (!response.ok) {
        const data = await response.json()
        throw new Error(data.detail || 'Failed to create job')
      }

      const data = await response.json()
      navigate(`/job/${data.job_id}`)
    } catch (err) {
      setError(err instanceof Error ? err.message : 'Unknown error')
    } finally {
      setSubmitting(false)
    }
  }

  return (
    <div className="max-w-2xl mx-auto">
      <div className="mb-8">
        <h1 className="text-3xl font-bold text-gray-900 dark:text-gray-100">New Annotation Job</h1>
        <p className="mt-2 text-gray-600 dark:text-gray-400">
          Upload your data and configure parameters for annotation
        </p>
      </div>

      {error && (
        <div className="mb-6 p-4 bg-red-50 dark:bg-red-900/30 border border-red-200 dark:border-red-800 rounded-lg flex items-center gap-2 text-red-700 dark:text-red-400">
          <AlertCircle className="w-5 h-5 flex-shrink-0" />
          {error}
        </div>
      )}

      <form onSubmit={handleSubmit} className="space-y-8">
        {/* File Uploads */}
        <div className="bg-white dark:bg-gray-900 rounded-lg border border-gray-200 dark:border-gray-700 p-6">
          <h2 className="text-lg font-semibold text-gray-900 dark:text-gray-100 mb-4 flex items-center gap-2">
            <Upload className="w-5 h-5" />
            Input Files
          </h2>

          <div className="space-y-4">
            {/* FASTA Upload */}
            <div>
              <label className="block text-sm font-medium text-gray-700 dark:text-gray-300 mb-2 flex items-center">
                Background Proteome FASTA *
                <Tooltip text="A FASTA file containing all protein sequences from the metaproteome community. Peptides are searched against this database." />
              </label>
              <div className="border-2 border-dashed border-gray-300 dark:border-gray-600 rounded-lg p-4 text-center hover:border-indigo-400 dark:hover:border-indigo-500 transition-colors">
                <input
                  type="file"
                  accept=".fasta,.fa,.faa"
                  onChange={(e) => setFastaFile(e.target.files?.[0] || null)}
                  className="hidden"
                  id="fasta-upload"
                />
                <label htmlFor="fasta-upload" className="cursor-pointer">
                  {fastaFile ? (
                    <div className="flex items-center justify-center gap-2 text-indigo-600 dark:text-indigo-400">
                      <FileText className="w-5 h-5" />
                      {fastaFile.name}
                    </div>
                  ) : (
                    <div className="text-gray-500 dark:text-gray-400">
                      <Upload className="w-8 h-8 mx-auto mb-2" />
                      <span>Click to upload FASTA file</span>
                    </div>
                  )}
                </label>
              </div>
            </div>

            {/* Peptide Files Upload */}
            <div>
              <label className="block text-sm font-medium text-gray-700 dark:text-gray-300 mb-2 flex items-center">
                Peptide Lists (CSV/TSV) *
                <Tooltip text="Two-column files (peptide sequence and count/abundance), one per sample. Peptides are matched exactly against the background proteome." />
              </label>
              <div className="border-2 border-dashed border-gray-300 dark:border-gray-600 rounded-lg p-4 text-center hover:border-indigo-400 dark:hover:border-indigo-500 transition-colors">
                <input
                  type="file"
                  accept=".csv,.tsv,.txt"
                  multiple
                  onChange={(e) => setPeptideFiles(Array.from(e.target.files || []))}
                  className="hidden"
                  id="peptide-upload"
                />
                <label htmlFor="peptide-upload" className="cursor-pointer">
                  {peptideFiles.length > 0 ? (
                    <div className="text-indigo-600 dark:text-indigo-400">
                      <FileText className="w-5 h-5 mx-auto mb-2" />
                      {peptideFiles.length} file(s) selected
                      <ul className="text-sm mt-2">
                        {peptideFiles.map((f, i) => (
                          <li key={i}>{f.name}</li>
                        ))}
                      </ul>
                    </div>
                  ) : (
                    <div className="text-gray-500 dark:text-gray-400">
                      <Upload className="w-8 h-8 mx-auto mb-2" />
                      <span>Click to upload peptide files</span>
                    </div>
                  )}
                </label>
              </div>
            </div>
          </div>
        </div>

        {/* Parameters */}
        <div className="bg-white dark:bg-gray-900 rounded-lg border border-gray-200 dark:border-gray-700 p-6">
          <h2 className="text-lg font-semibold text-gray-900 dark:text-gray-100 mb-4 flex items-center gap-2">
            <Settings className="w-5 h-5" />
            Parameters
          </h2>

          <div className="grid grid-cols-2 gap-4">
            <div>
              <label className="block text-sm font-medium text-gray-700 dark:text-gray-300 mb-1 flex items-center">
                Search Tool
                <Tooltip text="The homology search engine used to match peptides against the background proteome." />
              </label>
              <select
                value={searchTool}
                onChange={(e) => setSearchTool(e.target.value)}
                className="w-full px-3 py-2 border border-gray-300 dark:border-gray-600 rounded-lg bg-white dark:bg-gray-800 text-gray-900 dark:text-gray-100 focus:ring-2 focus:ring-indigo-500 focus:border-indigo-500"
              >
                <option value="diamond">DIAMOND{diamondVersion ? ` v${diamondVersion}` : ''}</option>
              </select>
            </div>

            <div>
              <label className="block text-sm font-medium text-gray-700 dark:text-gray-300 mb-1 flex items-center">
                Annotated Database
                <Tooltip text="The annotated reference database (e.g. UniProt) to search against using DIAMOND. Database files (.dmnd) are configured by the server administrator." />
              </label>
              <select
                value={dbChoice}
                onChange={(e) => setDbChoice(e.target.value)}
                className="w-full px-3 py-2 border border-gray-300 dark:border-gray-600 rounded-lg bg-white dark:bg-gray-800 text-gray-900 dark:text-gray-100 focus:ring-2 focus:ring-indigo-500 focus:border-indigo-500"
              >
                {databases.map((db) => (
                  <option key={db.path} value={db.path} title={db.description}>
                    {db.name}
                  </option>
                ))}
              </select>
            </div>

            <div>
              <label className="block text-sm font-medium text-gray-700 dark:text-gray-300 mb-1 flex items-center">
                Max E-value
                <Tooltip text="Maximum e-value threshold for DIAMOND homology hits. Also used as a DIAMOND pre-filter. Lower values are more stringent (e.g., 1e-10 keeps only high-confidence alignments)." />
              </label>
              <input
                type="text"
                value={maxEvalue}
                onChange={(e) => setMaxEvalue(e.target.value)}
                placeholder="e.g., 1e-10"
                className="w-full px-3 py-2 border border-gray-300 dark:border-gray-600 rounded-lg bg-white dark:bg-gray-800 text-gray-900 dark:text-gray-100 placeholder-gray-400 dark:placeholder-gray-500 focus:ring-2 focus:ring-indigo-500 focus:border-indigo-500"
              />
            </div>

            <div>
              <label className="block text-sm font-medium text-gray-700 dark:text-gray-300 mb-1 flex items-center">
                Min % Identity
                <Tooltip text="Minimum percent identity (0–100) for retaining DIAMOND hits. Applied after the search as a post-filter. Hits below this threshold are discarded before annotation." />
              </label>
              <input
                type="text"
                value={minPident}
                onChange={(e) => setMinPident(e.target.value)}
                placeholder="e.g., 80"
                className="w-full px-3 py-2 border border-gray-300 dark:border-gray-600 rounded-lg bg-white dark:bg-gray-800 text-gray-900 dark:text-gray-100 placeholder-gray-400 dark:placeholder-gray-500 focus:ring-2 focus:ring-indigo-500 focus:border-indigo-500"
              />
            </div>

            <div>
              <label className="block text-sm font-medium text-gray-700 dark:text-gray-300 mb-1 flex items-center">
                Top K Hits
                <Tooltip text="Number of top-scoring DIAMOND hits to keep per query protein, ranked by bitscore. Tie-aware: if multiple hits share the same bitscore at the Kth position, all tied hits are retained to avoid arbitrary bias in annotation." />
              </label>
              <input
                type="text"
                value={topK}
                onChange={(e) => setTopK(e.target.value)}
                placeholder="e.g., 1"
                className="w-full px-3 py-2 border border-gray-300 dark:border-gray-600 rounded-lg bg-white dark:bg-gray-800 text-gray-900 dark:text-gray-100 placeholder-gray-400 dark:placeholder-gray-500 focus:ring-2 focus:ring-indigo-500 focus:border-indigo-500"
              />
            </div>
          </div>
        </div>

        {/* Enrichment Analysis */}
        <div className="bg-white dark:bg-gray-900 rounded-lg border border-gray-200 dark:border-gray-700 p-6">
          <h2 className="text-lg font-semibold text-gray-900 dark:text-gray-100 mb-4 flex items-center gap-2">
            <FlaskConical className="w-5 h-5" />
            Enrichment Analysis
          </h2>

          <div className="space-y-4">
            <div className="flex items-start gap-3">
              <input
                type="checkbox"
                id="compute-enrichment"
                checked={computeEnrichment}
                onChange={(e) => setComputeEnrichment(e.target.checked)}
                className="mt-1 h-4 w-4 text-indigo-600 focus:ring-indigo-500 border-gray-300 dark:border-gray-600 rounded"
              />
              <div>
                <label htmlFor="compute-enrichment" className="text-sm font-medium text-gray-700 dark:text-gray-300 flex items-center cursor-pointer">
                  Calculate enrichment p-values
                  <Tooltip text="Uses Monte Carlo randomization to test whether GO terms are unusually enriched or depleted within specific taxa. Adds p-values and BH-adjusted q-values to the GO × taxonomy cross-tabulation report. This is compute-intensive and will increase job runtime." />
                </label>
                <p className="text-xs text-gray-500 dark:text-gray-400 mt-0.5">
                  Adds significance testing to GO × taxonomy combinations
                </p>
              </div>
            </div>

            {computeEnrichment && (
              <div className="ml-7">
                <label className="block text-sm font-medium text-gray-700 dark:text-gray-300 mb-1 flex items-center">
                  Monte Carlo Iterations
                  <Tooltip text="Number of randomization iterations. More iterations give more precise p-values but take longer. 1,000 is suitable for routine analysis; use 10,000 for high-resolution significance." />
                </label>
                <input
                  type="number"
                  min="100"
                  max="10000"
                  step="100"
                  value={enrichmentIterations}
                  onChange={(e) => setEnrichmentIterations(e.target.value)}
                  className="w-32 px-3 py-2 border border-gray-300 dark:border-gray-600 rounded-lg bg-white dark:bg-gray-800 text-gray-900 dark:text-gray-100 focus:ring-2 focus:ring-indigo-500 focus:border-indigo-500"
                />
                <p className="text-xs text-gray-500 dark:text-gray-400 mt-1">Range: 100 – 10,000</p>
              </div>
            )}
          </div>
        </div>

        {/* Notification */}
        <div className="bg-white dark:bg-gray-900 rounded-lg border border-gray-200 dark:border-gray-700 p-6">
          <h2 className="text-lg font-semibold text-gray-900 dark:text-gray-100 mb-4 flex items-center gap-2">
            Notification
          </h2>
          <div>
            <label className="block text-sm font-medium text-gray-700 dark:text-gray-300 mb-1 flex items-center">
              Email Address (optional)
              <Tooltip text="If provided, a notification email will be sent when your job completes or fails. The email will include your uploaded filenames, chosen parameters, and a link to view results." />
            </label>
            <input
              type="email"
              value={notificationEmail}
              onChange={(e) => setNotificationEmail(e.target.value)}
              placeholder="user@example.com"
              className="w-full px-3 py-2 border border-gray-300 dark:border-gray-600 rounded-lg bg-white dark:bg-gray-800 text-gray-900 dark:text-gray-100 placeholder-gray-400 dark:placeholder-gray-500 focus:ring-2 focus:ring-indigo-500 focus:border-indigo-500"
            />
          </div>
        </div>

        {/* Submit */}
        <div className="flex justify-end">
          <button
            type="submit"
            disabled={submitting}
            className="inline-flex items-center gap-2 px-6 py-3 bg-indigo-600 text-white rounded-lg hover:bg-indigo-700 disabled:opacity-50 disabled:cursor-not-allowed transition-colors"
          >
            {submitting ? (
              <>
                <Loader2 className="w-5 h-5 animate-spin" />
                Submitting...
              </>
            ) : (
              <>
                <Upload className="w-5 h-5" />
                Submit Job
              </>
            )}
          </button>
        </div>
      </form>
    </div>
  )
}
