import { Link } from 'react-router-dom'
import { Dna, Plus } from 'lucide-react'
import { useEffect, useState } from 'react'

interface LayoutProps {
  children: React.ReactNode
}

export default function Layout({ children }: LayoutProps) {
  const [version, setVersion] = useState<string>('')

  useEffect(() => {
    fetch('/api/version')
      .then(res => res.json())
      .then(data => setVersion(data.version))
      .catch(() => setVersion(''))
  }, [])

  return (
    <div className="min-h-screen flex flex-col">
      <header className="bg-white border-b border-gray-200 shadow-sm">
        <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8">
          <div className="flex justify-between items-center h-16">
            <Link to="/" className="flex items-center gap-2 text-xl font-bold text-indigo-600">
              <Dna className="w-8 h-8" />
              <span>Metagomics 2</span>
              {version && <span className="text-sm font-normal text-gray-500">v{version}</span>}
            </Link>
            <nav className="flex items-center gap-4">
              <Link
                to="/new"
                className="inline-flex items-center gap-2 px-4 py-2 bg-indigo-600 text-white rounded-lg hover:bg-indigo-700 transition-colors"
              >
                <Plus className="w-4 h-4" />
                New Job
              </Link>
            </nav>
          </div>
        </div>
      </header>
      <main className="flex-1">
        <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8 py-8">
          {children}
        </div>
      </main>
      <footer className="bg-white border-t border-gray-200 py-4">
        <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8 text-center text-gray-500 text-sm">
          Metagomics 2 - Metaproteomics Annotation Tool
          {version && <span className="ml-2">· Version {version}</span>}
        </div>
      </footer>
    </div>
  )
}
