import { describe, it, expect, beforeEach } from 'vitest'
import { render, screen, fireEvent } from '@testing-library/react'
import ThemeToggle from '../ThemeToggle'
import { ThemeProvider } from '../../ThemeContext'

function renderWithTheme() {
  return render(
    <ThemeProvider>
      <ThemeToggle />
    </ThemeProvider>
  )
}

describe('ThemeToggle', () => {
  beforeEach(() => {
    localStorage.clear()
    document.documentElement.classList.remove('dark')
  })

  it('renders a button with accessible label', () => {
    renderWithTheme()
    const button = screen.getByRole('button', { name: /switch to dark mode/i })
    expect(button).toBeTruthy()
  })

  it('toggles to dark mode on click', () => {
    renderWithTheme()
    const button = screen.getByRole('button', { name: /switch to dark mode/i })
    fireEvent.click(button)

    expect(document.documentElement.classList.contains('dark')).toBe(true)
    // After toggling to dark, the label should indicate switching to light
    expect(screen.getByRole('button', { name: /switch to light mode/i })).toBeTruthy()
  })

  it('toggles back to light mode on second click', () => {
    renderWithTheme()
    const button = screen.getByRole('button', { name: /switch to dark mode/i })

    fireEvent.click(button)
    expect(document.documentElement.classList.contains('dark')).toBe(true)

    fireEvent.click(screen.getByRole('button', { name: /switch to light mode/i }))
    expect(document.documentElement.classList.contains('dark')).toBe(false)
  })

  it('starts in dark mode if localStorage has dark theme', () => {
    localStorage.setItem('metagomics-theme', 'dark')
    renderWithTheme()
    expect(screen.getByRole('button', { name: /switch to light mode/i })).toBeTruthy()
    expect(document.documentElement.classList.contains('dark')).toBe(true)
  })
})
