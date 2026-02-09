import { Suspense, lazy } from 'react'
import { Routes, Route } from 'react-router-dom'
import Layout from './components/Layout'
import NewJobPage from './pages/NewJobPage'
import JobPage from './pages/JobPage'
import AdminPage from './pages/AdminPage'
import GoDagPage from './pages/GoDagPage'

const TaxonomyPage = lazy(() => import('./pages/TaxonomyPage'))

function App() {
  return (
    <Layout>
      <Suspense fallback={<div className="flex items-center justify-center py-12 text-gray-500">Loading…</div>}>
        <Routes>
          <Route path="/" element={<NewJobPage />} />
          <Route path="/job/:jobId" element={<JobPage />} />
          <Route path="/job/:jobId/go/:listId" element={<GoDagPage />} />
          <Route path="/job/:jobId/taxonomy/:listId" element={<TaxonomyPage />} />
          <Route path="/admin" element={<AdminPage />} />
        </Routes>
      </Suspense>
    </Layout>
  )
}

export default App
