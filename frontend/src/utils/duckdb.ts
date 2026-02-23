import * as duckdb from '@duckdb/duckdb-wasm'

interface DuckDBInstance {
  db: duckdb.AsyncDuckDB
  conn: duckdb.AsyncDuckDBConnection
}

let instance: DuckDBInstance | null = null
let initPromise: Promise<DuckDBInstance> | null = null

export async function getDuckDB(): Promise<DuckDBInstance> {
  if (instance) return instance
  if (initPromise) return initPromise

  initPromise = (async () => {
    const JSDELIVR_BUNDLES = duckdb.getJsDelivrBundles()
    const bundle = await duckdb.selectBundle(JSDELIVR_BUNDLES)

    const workerUrl = URL.createObjectURL(
      new Blob([`importScripts("${bundle.mainWorker!}");`], { type: 'text/javascript' })
    )
    const worker = new Worker(workerUrl)
    const logger = new duckdb.ConsoleLogger()
    const db = new duckdb.AsyncDuckDB(logger, worker)
    await db.instantiate(bundle.mainModule, bundle.pthreadWorker)
    URL.revokeObjectURL(workerUrl)

    const conn = await db.connect()
    instance = { db, conn }
    return instance
  })()

  return initPromise
}

export async function registerMappingFile(jobId: string, listId: string): Promise<void> {
  const { db, conn } = await getDuckDB()
  const url = `${window.location.origin}/api/jobs/${jobId}/results/${listId}/peptide_mapping.parquet`
  await db.registerFileURL('mappings.parquet', url, duckdb.DuckDBDataProtocol.HTTP, false)
  await conn.query(`
    CREATE OR REPLACE VIEW mappings AS
    SELECT * FROM read_parquet('mappings.parquet')
  `)
}
