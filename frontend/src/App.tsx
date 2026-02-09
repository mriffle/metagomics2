import { Routes, Route } from 'react-router-dom'
import Layout from './components/Layout'
import NewJobPage from './pages/NewJobPage'
import JobPage from './pages/JobPage'
import AdminPage from './pages/AdminPage'
import GoDagPage from './pages/GoDagPage'

function App() {
  return (
    <Layout>
      <Routes>
        <Route path="/" element={<NewJobPage />} />
        <Route path="/job/:jobId" element={<JobPage />} />
        <Route path="/job/:jobId/go/:listId" element={<GoDagPage />} />
        <Route path="/admin" element={<AdminPage />} />
      </Routes>
    </Layout>
  )
}

export default App
