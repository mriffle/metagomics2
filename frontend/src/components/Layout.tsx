import { Link } from 'react-router-dom'
import { useEffect, useState } from 'react'
import MetagomicsLogo from './MetagomicsLogo'
import ThemeToggle from './ThemeToggle'

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
    <div className="min-h-screen flex flex-col bg-gray-50 dark:bg-gray-950 text-gray-900 dark:text-gray-100">
      <header className="bg-white dark:bg-gray-900 border-b border-gray-200 dark:border-gray-700 shadow-sm">
        <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8">
          <div className="flex justify-between items-center h-16">
            <Link to="/" className="flex items-center gap-2 text-xl font-bold text-indigo-600 dark:text-indigo-400">
              <MetagomicsLogo className="h-10 w-auto" />
              {version && <span className="text-sm font-normal text-gray-500 dark:text-gray-400">v{version}</span>}
            </Link>
            <nav className="flex items-center gap-4">
              <ThemeToggle />
              <Link
                to="/admin"
                className="text-sm text-gray-500 hover:text-gray-700 dark:text-gray-400 dark:hover:text-gray-200 transition-colors"
              >
                Admin
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
      <footer className="bg-white dark:bg-gray-900 border-t border-gray-200 dark:border-gray-700 py-4">
        <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8 text-center text-gray-500 dark:text-gray-400 text-sm">
          Metagomics 2 - Metaproteomics Annotation Tool
          {version && <span className="ml-2">· Version {version}</span>}
        </div>
      </footer>
    </div>
  )
}
