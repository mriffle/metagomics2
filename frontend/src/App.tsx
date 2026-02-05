import { Routes, Route } from 'react-router-dom'
import Layout from './components/Layout'
import HomePage from './pages/HomePage'
import NewJobPage from './pages/NewJobPage'
import JobPage from './pages/JobPage'

function App() {
  return (
    <Layout>
      <Routes>
        <Route path="/" element={<HomePage />} />
        <Route path="/new" element={<NewJobPage />} />
        <Route path="/job/:jobId" element={<JobPage />} />
      </Routes>
    </Layout>
  )
}

export default App
