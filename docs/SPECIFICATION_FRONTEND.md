# Metagomics 2 — Frontend Specification

## 1. Purpose

The Metagomics 2 frontend is a single-page application (SPA) that provides a web interface for submitting metaproteomics annotation jobs and interactively visualizing results. Users can:

1. **Submit jobs**: Upload a background proteome FASTA and peptide list files, configure analysis parameters, and submit to the backend pipeline
2. **Monitor jobs**: Track job progress in real time with auto-polling
3. **Download results**: Download individual CSV/Parquet result files or a ZIP of all results
4. **Visualize taxonomy**: Interactive sunburst, treemap, icicle, and Sankey charts of taxonomic assignments
5. **Visualize Gene Ontology**: Interactive DAG (directed acyclic graph) of GO term assignments using Cytoscape.js
6. **Cross-filter**: Filter taxonomy by GO term or GO by taxon using the combo cross-tabulation data
7. **Drill into peptide details**: Click a chart node to see which peptides, background proteins, and annotated proteins contribute, queried in real time from Parquet files via DuckDB-WASM
8. **Administer**: Password-protected admin dashboard to view all jobs

---

## 2. Developer Directives

These directives **must** be followed for all frontend development:

- **Do NOT install npm packages on the host**. All npm dependencies are managed via Docker. The `Dockerfile` `frontend-builder` stage handles `npm install` and `npm run build`.
- **TypeScript type safety**: All new frontend code must pass `tsc --noEmit`. Strict mode is enabled (`strict: true` in `tsconfig.json`).
- **Frontend tests run in Docker**:
  ```bash
  docker build --target frontend-builder -t metagomics2-frontend-test .
  docker run --rm metagomics2-frontend-test npx vitest run
  ```
- **TailwindCSS for styling**: Use Tailwind utility classes with `dark:` variants for dark mode. No custom CSS files beyond `index.css` (which imports Tailwind directives and defines CSS custom properties for theme tokens).
- **Dark mode support is mandatory**: All new UI elements must include `dark:` Tailwind variants. Visualization colors consumed by Cytoscape/Plotly JavaScript APIs must be defined in `utils/colors.ts` (the centralized color constants module), not inline. See Section 11 for the full dark mode architecture.
- **Lucide React for icons**: Use `lucide-react` for all icons. Do not add other icon libraries.
- **Write tests for new utility functions**: All parsers and data transformation functions in `utils/` have corresponding test files. Maintain this pattern.

---

## 3. Technology Stack

| Technology | Version | Purpose |
|------------|---------|---------|
| **React** | ^18.2 | UI framework |
| **React Router DOM** | ^6.21 | Client-side routing |
| **TypeScript** | ^5.2 | Type safety |
| **Vite** | ^5.0 | Build tool and dev server |
| **TailwindCSS** | ^3.4 | Utility-first CSS |
| **Plotly.js** | ^2.35 | Taxonomy charts (sunburst, treemap, icicle, Sankey) |
| **react-plotly.js** | ^2.6 | React wrapper for Plotly (via factory pattern) |
| **Cytoscape.js** | ^3.30 | GO DAG graph visualization |
| **cytoscape-dagre** | ^2.5 | Hierarchical DAG layout for Cytoscape |
| **cytoscape-svg** | ^0.4 | SVG export for Cytoscape graphs |
| **DuckDB-WASM** | ^1.29 | In-browser SQL engine for querying Parquet files |
| **Lucide React** | ^0.303 | Icon library |
| **Vitest** | ^2.0 | Test runner |
| **Testing Library** | ^16.0 (react), ^6.4 (jest-dom), ^14.5 (user-event) | Component testing |
| **jsdom** | ^24.0 | DOM environment for tests |
| **PostCSS** | ^8.4 | CSS processing (Tailwind plugin) |
| **Autoprefixer** | ^10.4 | Vendor prefixing |

---

## 4. Project Structure

```
frontend/
├── index.html                         # HTML entry point (mounts #root, loads main.tsx)
├── package.json                       # Dependencies and scripts
├── vite.config.ts                     # Vite + Vitest config (proxy, test env)
├── tsconfig.json                      # TypeScript config (strict, ES2020, react-jsx)
├── tsconfig.node.json                 # TypeScript config for Vite config file
├── tailwind.config.js                 # Tailwind content paths + darkMode: 'class'
├── postcss.config.js                  # PostCSS plugins (tailwindcss, autoprefixer)
└── src/
    ├── main.tsx                       # React entry: ThemeProvider + BrowserRouter + App
    ├── App.tsx                        # Route definitions + Layout wrapper
    ├── ThemeContext.tsx                # Theme provider, context, and useTheme hook (light/dark)
    ├── index.css                      # Tailwind directives + CSS custom properties for theme tokens
    ├── test-setup.ts                  # Vitest setup (imports jest-dom matchers)
    ├── plotly.d.ts                    # Type declarations for plotly.js-dist-min, react-plotly.js/factory
    ├── cytoscape-dagre.d.ts           # Type declaration for cytoscape-dagre
    ├── cytoscape-svg.d.ts             # Type declaration for cytoscape-svg
    ├── __tests__/
    │   └── ThemeContext.test.tsx       # Tests for ThemeContext/ThemeProvider
    ├── pages/
    │   ├── NewJobPage.tsx             # Job submission form
    │   ├── JobPage.tsx                # Job status, progress, results, downloads
    │   ├── GoDagPage.tsx              # GO DAG visualization page
    │   ├── TaxonomyPage.tsx           # Taxonomy chart visualization page
    │   ├── AdminPage.tsx              # Admin login + job list dashboard
    │   └── HomePage.tsx               # Recent jobs list (currently unused in routing)
    ├── components/
    │   ├── Layout.tsx                 # App shell: header, nav, footer, theme toggle, version display
    │   ├── ThemeToggle.tsx            # Light/dark mode toggle switch (Sun/Moon icons)
    │   ├── GoDagViewer.tsx            # Cytoscape-based GO DAG renderer
    │   ├── GoDagControls.tsx          # GO DAG controls (namespace, metric, filter, export)
    │   ├── TaxonomyChart.tsx          # Plotly-based taxonomy chart renderer
    │   ├── TaxonomyControls.tsx       # Taxonomy controls (chart type, rank depth, filter, export)
    │   ├── PeptideDetailsPane.tsx     # DuckDB-powered peptide drill-down panel
    │   ├── Autocomplete.tsx           # Reusable autocomplete/search input
    │   ├── UniprotProteinLabel.tsx    # UniProt accession link with hover tooltip
    │   └── __tests__/
    │       ├── PeptideDetailsPane.test.tsx
    │       └── ThemeToggle.test.tsx    # Tests for ThemeToggle component
    └── utils/
        ├── csvParser.ts               # RFC-compliant CSV line parser
        ├── taxonomyParser.ts          # Taxonomy CSV parsing, canonical rank filtering, placeholder insertion
        ├── goParser.ts                # GO terms CSV parsing
        ├── comboParser.ts             # GO-taxonomy combo CSV parsing and reshaping
        ├── colors.ts                  # Centralized light/dark mode color constants for visualizations and badges
        ├── duckdb.ts                  # DuckDB-WASM singleton initialization and Parquet registration
        ├── uniprot.ts                 # UniProt accession extraction, URL building, REST API info fetching
        └── __tests__/
            ├── csvParser.test.ts
            ├── taxonomyParser.test.ts
            ├── goParser.test.ts
            └── comboParser.test.ts
```

---

## 5. Routing

Defined in `App.tsx`. All routes are wrapped in a `<Layout>` shell and `<Suspense>` for lazy-loaded pages.

| Route | Component | Description |
|-------|-----------|-------------|
| `/` | `NewJobPage` | Job submission form (home page) |
| `/job/:jobId` | `JobPage` | Job status, progress, results |
| `/job/:jobId/go/:listId` | `GoDagPage` | GO DAG visualization for a specific peptide list |
| `/job/:jobId/taxonomy/:listId` | `TaxonomyPage` | Taxonomy chart for a specific peptide list (lazy-loaded) |
| `/admin` | `AdminPage` | Admin dashboard (password-protected) |

**Note**: `TaxonomyPage` is lazy-loaded via `React.lazy()` because it imports the large Plotly.js library.

**Note**: `HomePage.tsx` exists in the codebase but is **not currently routed** — the root `/` path maps to `NewJobPage`.

---

## 6. Pages

### 6.1 NewJobPage (`pages/NewJobPage.tsx`)

The job submission form. Fetches server configuration on mount (`GET /api/config`) to populate the database dropdown and DIAMOND version display.

**State**:
- `fastaFile: File | null` — Selected FASTA file
- `peptideFiles: File[]` — Selected peptide list files (multiple)
- `searchTool: string` — Always `"diamond"`
- `dbChoice: string` — Selected annotated database path
- `maxEvalue, minPident, topK: string` — Filter parameters (text inputs)
- `notificationEmail: string` — Optional email for notifications
- `computeEnrichment: boolean` — Whether to compute Monte Carlo enrichment p-values
- `enrichmentIterations: string` — Number of Monte Carlo iterations (100–10000, default 1000)

**Submission flow**:
1. Validates FASTA and peptide files are selected
2. Builds `FormData` with `fasta`, `peptides` (repeatable), and `params` (JSON string)
3. `POST /api/jobs` with multipart body
4. On success: `navigate(/job/${data.job_id})`

**UI features**:
- Drag-style file upload areas with visual feedback
- Inline `Tooltip` component for parameter help text
- Error display with `AlertCircle` icon
- Loading spinner during submission

### 6.2 JobPage (`pages/JobPage.tsx`)

Displays job status, parameters, progress, and results.

**Polling**: Fetches `GET /api/jobs/:jobId` every 3 seconds. Stops polling when status is `completed` or `failed`.

**Sections**:
- **Header**: Status icon, job ID (monospace), "Change Hash" button, creation time
- **Parameters**: Inline display of FASTA filename, database, filter params
- **Progress bar**: Shown for `running`/`queued` jobs. Shows `current_step` (left) and overall completion percentage (right). Progress uses a weighted milestone system (0–1000 scale) that reflects all pipeline stages — initialization, peptide matching, homology search, filtering, and per-list processing — not just peptide list completion. See the backend spec's `PipelineProgress` section for milestone weights.
- **Error message**: Shown for `failed` jobs
- **Peptide Lists**: Each list shows filename, peptide count, match count, status badge
- **Per-list results** (when `completed`):
  - Download links for each result file (with tooltips explaining each file)
  - "View results" links to GO DAG and Taxonomy visualization pages
  - "Download All Results (ZIP)" button

**Job ID regeneration**: "Change Hash" button calls `POST /api/jobs/:jobId/regenerate-id`, then navigates to the new URL with `replace: true`.

### 6.3 GoDagPage (`pages/GoDagPage.tsx`)

Interactive Gene Ontology DAG visualization.

**Data flow**:
1. Fetches `go_terms.csv` → `parseGoTermsCsv()` → `GoTermNode[]`
2. Fetches `taxonomy_nodes.csv` → `parseTaxonomyCsv()` → builds taxon autocomplete options
3. Optionally fetches `go_taxonomy_combo.csv` when a taxon filter is selected → `parseComboCsv()` → `comboRowsToGoTermNodes()`

**Filtering pipeline**:
1. Filter by namespace (`biological_process` | `cellular_component` | `molecular_function`)
2. Filter by `minRatioTotal` cutoff
3. Optionally filter by taxon (uses combo data)

**Child components**:
- `GoDagControls` — namespace tabs, metric selector, color picker, abundance cutoff, taxon autocomplete, export buttons
- `GoDagViewer` — Cytoscape.js graph renderer
- `PeptideDetailsPane` — peptide drill-down (shown when a GO node is clicked)

**Metric options**: `quantity`, `ratioTotal`, `ratioAnnotated`, `nPeptides`, plus `fractionOfTaxon`, `fractionOfGo`, and `qvalueGoForTaxon` when a taxon filter is active. The `qvalueGoForTaxon` metric uses inverted `-log10(q)` normalization so that low q-values (significant) appear as intense colors.

### 6.4 TaxonomyPage (`pages/TaxonomyPage.tsx`)

Interactive taxonomy visualization with multiple chart types.

**Data flow**:
1. Fetches `taxonomy_nodes.csv` → `parseTaxonomyCsv()` → `TaxonNode[]`
2. Runs through the **taxonomy processing pipeline**:
   - `filterCanonicalRanks()` — Keep only canonical ranks, re-link parents
   - `validateCanonicalHierarchy()` — Validate parent-child rank ordering
   - `filterByMaxRank()` — Limit depth to selected rank
   - `ensureStrictRankLayers()` — Insert placeholder nodes for missing intermediate ranks (except for Sankey)
3. Fetches `go_terms.csv` → builds GO autocomplete options
4. Optionally fetches `go_taxonomy_combo.csv` when a GO filter is selected → `comboRowsToTaxonNodes()`

**Chart types**: `sunburst`, `treemap`, `icicle`, `sankey`

**Child components**:
- `TaxonomyControls` — chart type buttons, rank depth selector, GO term autocomplete, abundance cutoff, export buttons
- `TaxonomyChart` — Plotly.js chart renderer
- `PeptideDetailsPane` — peptide drill-down (shown when a taxonomy node is clicked)

**Node click behavior**: Clicking a taxonomy node selects it and all its descendants (via `getDescendantTaxIds()`), which is passed to `PeptideDetailsPane` for filtering.

### 6.5 AdminPage (`pages/AdminPage.tsx`)

Password-protected admin dashboard.

**Two sub-components**:
- `LoginForm` — Password input, calls `POST /api/admin/auth`, stores token in `sessionStorage`
- `JobList` — Polls `GET /api/admin/jobs` every 5 seconds with `Authorization: Bearer <token>`. Shows a table of all jobs with links, timestamps, status badges, and progress.

Auto-logs out if the server returns 401.

### 6.6 HomePage (`pages/HomePage.tsx`)

Lists recent jobs (polls `GET /api/jobs` every 5 seconds). Contains a custom inline SVG Dna icon. **Not currently routed** — exists but unused.

---

## 7. Components

### 7.1 Layout (`components/Layout.tsx`)

App shell wrapping all pages. Fetches version from `GET /api/version` on mount.

- **Header**: Logo (Lucide `Dna` icon) + "Metagomics 2" link to `/`, version badge, `ThemeToggle` switch, "Admin" nav link
- **Main**: `max-w-7xl` centered content area
- **Footer**: App name + version

All Layout elements use `dark:` Tailwind variants for dark mode styling (e.g., `bg-white dark:bg-gray-900`).

### 7.1a ThemeToggle (`components/ThemeToggle.tsx`)

Pill-shaped toggle switch for switching between light and dark mode. Rendered in the Layout header between the version badge and the Admin link.

- Uses Lucide `Sun` and `Moon` icons
- Sliding knob animates left/right via CSS `translate-x` transition
- Background icons indicate the inactive mode at reduced opacity
- Accessible `aria-label` updates based on current state ("Switch to dark mode" / "Switch to light mode")
- Calls `toggleTheme()` from `useTheme()` context hook on click

### 7.2 GoDagViewer (`components/GoDagViewer.tsx`)

Renders the GO DAG using Cytoscape.js with the `dagre` layout.

**Key behaviors**:
- **Layout**: Top-to-bottom DAG using `dagre` layout engine with `tight-tree` ranker
- **Node coloring**: Based on selected metric, normalized to [0, 1] via min-max scaling. Log transform for `quantity` and `nPeptides`. 5-stop color ramp whose endpoints are determined by the current theme (defined in `GO_DAG` constants from `utils/colors.ts`). In dark mode, nodes have a subtle indigo glow via Cytoscape `shadow-*` styles. Node text color is chosen by computing the WCAG relative luminance of the actual interpolated background color and selecting light or dark text for contrast (defined as `textLight`/`textDark` in the `GO_DAG` constants).
- **Hover**: BFS traversal highlights the hovered node and **all ancestor nodes and edges up to the root** (not just immediate neighbors), using `cy.incomers('edge')` recursively. Shows a fixed-position tooltip with node details (ID, name, quantity, ratios, fractions). Tooltip styling adapts to dark mode.
- **Tooltip positioning**: Uses `useLayoutEffect` to reposition within viewport bounds, flipping horizontally/vertically if overflowing.
- **Click**: Notifies parent via `onNodeClick(nodeId)`. Background click clears selection.
- **Export**: PNG via `cy.png()`, SVG via `cy.svg()`. Functions are attached to the container DOM element for parent access (`__exportPng`, `__exportSvg`).
- **Metric update**: When metric or baseColor changes, colors are updated in a `cy.batch()` without re-running layout.

### 7.3 GoDagControls (`components/GoDagControls.tsx`)

Controls for the GO DAG visualization.

**Sections**:
- **Namespace tabs**: Toggle between GO namespaces
- **Metric selector**: Dropdown to choose coloring metric
- **Color picker**: Click the gradient bar to open a hidden `<input type="color">`
- **Export buttons**: PNG and SVG
- **Taxonomy filter**: `Autocomplete` component to filter GO terms by a specific taxon
- **Abundance cutoff**: Preset buttons (None, 0.01%, 0.1%, 1%, 5%, 10%) + custom percentage input

### 7.4 TaxonomyChart (`components/TaxonomyChart.tsx`)

Renders taxonomy data using Plotly.js via `react-plotly.js/factory`.

**Chart types**:
- **Sunburst**: Hierarchical pie chart. `branchvalues='total'`. Node colors assigned by domain ancestry with depth-based lightening.
- **Treemap**: Hierarchical rectangles with `squarify` packing and pathbar navigation.
- **Icicle**: Vertical partition chart with pathbar.
- **Sankey**: Flow diagram with explicit x-position columns for each rank. Semi-transparent link colors. No click/hover handlers (Sankey has built-in Plotly hover).

**Coloring strategy**: Each domain (top-level taxon below root) gets a distinct color from the domain palette defined in `utils/colors.ts`. In light mode, `DOMAIN_PALETTE_LIGHT` provides medium-saturation colorblind-friendly colors that are progressively lightened at depth. In dark mode, `DOMAIN_PALETTE_DARK` provides vivid high-saturation neon colors with minimal depth-based washout, and `textfont.color` is set to dark (`#111827`) for contrast against bright backgrounds. Root is white in light mode, dark gray (`rgb(17,24,39)`) in dark mode. All chart-specific color values (root, fallback, line separators, text, pathbar text) are read from `TAX_CHART` and `SANKEY` constants in `utils/colors.ts`.

**Hover tooltip**: Custom fixed-position tooltip (same pattern as GoDagViewer) showing tax ID, name, rank, quantity, ratios, fractions.

**Click**: Notifies parent via `onNodeClick(taxId)`. Not active for Sankey.

### 7.5 TaxonomyControls (`components/TaxonomyControls.tsx`)

Controls for taxonomy visualization.

**Sections**:
- **Chart type toggle**: Sunburst, Treemap, Icicle, Sankey buttons
- **Rank depth selector**: Dropdown from Domain to Species
- **Export buttons**: PNG and SVG (uses Plotly's `downloadImage`)
- **GO term filter**: `Autocomplete` to filter taxonomy by a specific GO term
- **Abundance cutoff**: Same preset + custom pattern as GoDagControls

### 7.6 PeptideDetailsPane (`components/PeptideDetailsPane.tsx`)

Drill-down panel showing peptide-level detail for a selected taxonomy node or GO term.

**Data engine**: Uses DuckDB-WASM to query `peptide_mapping.parquet` directly in the browser.

**Query flow**:
1. On mount (or when `jobId`/`listId` change): Registers the Parquet file URL as a DuckDB view via `registerMappingFile()`
2. On selection change: Builds a SQL `WHERE` clause:
   - `selectedTaxIds` → `list_has_any(peptide_lca_tax_ids, [ids...])`
   - `selectedGoId` → `list_contains(peptide_go_terms, 'GO:...')`
3. Executes the query and builds a three-level hierarchy: `peptide → background_protein → annotated_protein`

**UI hierarchy** (expandable tree):
- **Peptide** (monospace, with count of background proteins)
  - **Background protein** (indigo, with count of subjects)
    - **Annotated protein** (green, with e-value and % identity)
    - If the annotated protein ID is a UniProt accession, renders `UniprotProteinLabel` with hover tooltip

**CSV download**: Generates a CSV of all visible peptide details.

### 7.7 Autocomplete (`components/Autocomplete.tsx`)

Reusable search/autocomplete dropdown.

**Features**:
- Case-insensitive substring filtering
- Keyboard navigation (Arrow Up/Down, Enter, Escape)
- Max 50 results shown
- Clear button (X icon)
- Auto-scrolls highlighted item into view
- Closes on outside click
- Syncs input text with external `value` prop

### 7.8 UniprotProteinLabel (`components/UniprotProteinLabel.tsx`)

Renders a UniProt accession as a clickable link to uniprot.org with a hover tooltip.

**Behavior**:
1. Extracts accession from raw ID (handles `sp|ACC|NAME` and bare formats)
2. Fetches protein info from `https://rest.uniprot.org/uniprotkb/{ACC}?format=json` (8s timeout)
3. Caches results in a module-level `Map` (with in-flight deduplication)
4. On hover: Shows tooltip with reviewed/unreviewed status (Swiss-Prot vs TrEMBL), full protein name, gene name, organism
5. Link opens UniProt page in new tab

---

## 8. Utility Modules

### 8.1 csvParser (`utils/csvParser.ts`)

Single function `parseCSVLine(line: string): string[]` — parses one CSV line handling:
- Quoted fields (double-quote delimited)
- Escaped quotes (`""` → `"`)
- Commas inside quoted fields
- Mixed quoted/unquoted fields

Used by all other parsers to handle fields that may contain commas (e.g., species names like "Homo sapiens, neanderthalensis").

### 8.2 taxonomyParser (`utils/taxonomyParser.ts`)

**Types**:
- `TaxonNode` — `{ taxId, name, rank, parentTaxId, quantity, ratioTotal, ratioAnnotated, nPeptides, fractionOfTaxon?, fractionOfGo?, qvalueTaxonForGo?, qvalueGoForTaxon? }`
- `CanonicalRank` — `'root' | 'domain' | 'kingdom' | 'phylum' | 'class' | 'order' | 'family' | 'genus' | 'species'`

**Constants**:
- `CANONICAL_RANKS_ORDERED` — Array of 9 canonical ranks from root to species
- `CANONICAL_RANKS` — `Set<string>` for O(1) membership tests

**Functions**:

| Function | Purpose |
|----------|---------|
| `parseTaxonomyCsv(text)` | Parse `taxonomy_nodes.csv` into `TaxonNode[]` |
| `filterCanonicalRanks(nodes)` | Remove non-canonical ranks; re-link parents through non-canonical intermediates; normalize NCBI root ("no rank") to "root" |
| `validateCanonicalHierarchy(nodes)` | Check that every non-root node's parent exists and is at a strictly higher rank. Returns error messages. Rank gaps are allowed. |
| `filterByMaxRank(nodes, maxRank)` | Keep only nodes at or above the specified rank depth |
| `ensureStrictRankLayers(nodes)` | Insert placeholder nodes for missing intermediate ranks so chart ring depth matches rank depth. Placeholders have IDs like `__placeholder_kingdom_<childTaxId>` and names like `(no kingdom for Haptophyta)`. |
| `getDescendantTaxIds(rootTaxId, allNodes)` | BFS to collect a node and all its descendants. Used for peptide detail filtering. |

**Taxonomy processing pipeline** (used by TaxonomyPage):
```
parseTaxonomyCsv → filterCanonicalRanks → validateCanonicalHierarchy → filterByMaxRank → ensureStrictRankLayers → chart
```

**Key design decisions**:
- NCBI root node (tax_id 1) has rank "no rank" — `filterCanonicalRanks` normalizes this to "root"
- Not all NCBI lineages have every canonical rank (e.g., Bacteria has no kingdom) — `validateCanonicalHierarchy` allows rank gaps
- `ensureStrictRankLayers` inserts placeholder nodes so Plotly sunburst rings correspond to consistent rank levels
- Sankey charts skip `ensureStrictRankLayers` because they use explicit x-position columns for rank alignment

### 8.3 goParser (`utils/goParser.ts`)

**Types**:
- `GoTermNode` — `{ id, name, namespace, parentIds: string[], quantity, ratioTotal, ratioAnnotated, nPeptides, fractionOfTaxon?, fractionOfGo?, qvalueGoForTaxon?, qvalueTaxonForGo? }`

**Functions**:
- `parseGoTermsCsv(text)` — Parse `go_terms.csv` into `GoTermNode[]`. Parent IDs are semicolon-delimited.

### 8.4 comboParser (`utils/comboParser.ts`)

Handles the `go_taxonomy_combo.csv` cross-tabulation data for cross-filtering.

**Types**:
- `ComboRow` — `{ taxId, taxName, taxRank, parentTaxId, goId, goName, goNamespace, parentGoIds, quantity, fractionOfTaxon, fractionOfGo, ratioTotalTaxon, ratioTotalGo, nPeptides, pvalueGoForTaxon?, pvalueTaxonForGo?, qvalueGoForTaxon?, qvalueTaxonForGo? }`

**Functions**:
| Function | Purpose |
|----------|---------|
| `parseComboCsv(text)` | Parse combo CSV into `ComboRow[]` |
| `comboRowsToTaxonNodes(rows, goId)` | Filter by GO ID, reshape into `TaxonNode[]`. Uses `ratioTotalTaxon` as `ratioTotal`. |
| `comboRowsToGoTermNodes(rows, taxId)` | Filter by tax ID, reshape into `GoTermNode[]`. Uses `ratioTotalGo` as `ratioTotal`. |

### 8.5 duckdb (`utils/duckdb.ts`)

Singleton DuckDB-WASM instance for querying Parquet files in the browser.

**Functions**:
- `getDuckDB()` — Lazy-initializes a DuckDB-WASM instance using jsDelivr bundles. Returns `{ db, conn }`. The instance is created once and reused.
- `registerMappingFile(jobId, listId)` — Registers the `peptide_mapping.parquet` file for a specific job/list as an HTTP-backed Parquet source. Creates a `mappings` view for SQL queries.

**Architecture**: DuckDB-WASM runs entirely in the browser via WebAssembly. It fetches Parquet files via HTTP range requests, enabling efficient column-oriented queries without downloading the entire file.

### 8.6 uniprot (`utils/uniprot.ts`)

UniProt protein annotation utilities.

**Functions**:
| Function | Purpose |
|----------|---------|
| `extractUniprotAccession(id)` | Extract accession from `sp\|ACC\|NAME`, `tr\|ACC\|NAME`, or bare format. Returns `null` if not a valid UniProt accession. |
| `isUniprotAccession(id)` | Boolean convenience wrapper |
| `uniprotUrl(accession)` | Build UniProt KB URL |
| `fetchUniprotInfo(accession)` | Fetch protein metadata from UniProt REST API. Returns `UniprotInfo` or `null`. Module-level cache + in-flight deduplication. 8-second timeout. |

**UniprotInfo type**: `{ accession, name, fullName, organism, gene, reviewed }`

**Accession regex**: Matches standard UniProt formats: `[OPQ][0-9][A-Z0-9]{3}[0-9]` or `[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}` with optional isoform suffix (`-\d+`).

### 8.7 colors (`utils/colors.ts`)

Centralized color definitions for light and dark mode. All visualization colors consumed by Cytoscape.js and Plotly.js APIs are defined here as TypeScript constants rather than in CSS, because these libraries require colors passed as JavaScript values to their programmatic APIs (they render to `<canvas>`/SVG, not DOM elements styled by CSS classes).

**Exported constants**:

| Constant | Purpose |
|----------|---------|
| `DOMAIN_PALETTE_LIGHT` | Light-mode taxonomy domain palette (10 colorblind-friendly colors as `number[][]`) |
| `DOMAIN_PALETTE_DARK` | Dark-mode taxonomy domain palette (10 vivid neon colors as `number[][]`) |
| `TAX_CHART` | Taxonomy chart colors keyed by `light`/`dark` (root color, fallback color, line separator color, inside-text color, pathbar text color) |
| `GO_DAG` | GO DAG Cytoscape colors keyed by `light`/`dark` (node border width/color, background, text color, shadow/glow parameters, edge color, highlight border/edge/shadow, color ramp blend parameters, light/dark text colors for contrast, export background) |
| `SANKEY` | Sankey-specific node line colors keyed by `light`/`dark` |
| `PLOTLY_LAYOUT` | Plotly layout font colors keyed by `light`/`dark` |
| `STATUS_BADGE_CLASSES` | `Record<string, string>` mapping status names (`uploaded`, `queued`, `pending`, `running`, `completed`, `done`, `failed`) to Tailwind class strings with both light and `dark:` variants |
| `STATUS_BADGE_DEFAULT` | Fallback badge class string for unknown status values |

**Design rationale**: Cytoscape and Plotly render to `<canvas>` or their own SVG, bypassing DOM CSS styling entirely. Colors must be passed as JS values to their APIs. Additionally, taxonomy domain colors require runtime blending (computed per-node based on domain ancestry and depth), which cannot be expressed as static CSS classes. Status badge classes are Tailwind strings composed dynamically in JSX template literals, so they are also stored as JS constants. Keeping all of these in one module ensures a single source of truth for the dark mode color palette.

---

## 9. API Communication

The frontend communicates with the backend exclusively through the `/api/` prefix. In development, Vite proxies `/api` to `http://localhost:8000`. In production, the FastAPI server serves both the API and the built SPA.

### API Endpoints Used by Frontend

| Endpoint | Method | Used By | Purpose |
|----------|--------|---------|---------|
| `/api/version` | GET | `Layout` | Fetch version string |
| `/api/config` | GET | `NewJobPage` | Fetch available databases, DIAMOND version |
| `/api/jobs` | POST | `NewJobPage` | Submit new job (multipart form) |
| `/api/jobs/:id` | GET | `JobPage`, `GoDagPage`, `TaxonomyPage` | Get job status and peptide list info |
| `/api/jobs/:id/regenerate-id` | POST | `JobPage` | Generate new job URL |
| `/api/jobs/:id/results/:listId/:filename` | GET | `JobPage`, `GoDagPage`, `TaxonomyPage`, `PeptideDetailsPane` | Download/fetch result files |
| `/api/jobs/:id/results/all_results.zip` | GET | `JobPage` | Download all results as ZIP |
| `/api/admin/auth` | POST | `AdminPage` | Admin login |
| `/api/admin/jobs` | GET | `AdminPage` | List all jobs (requires Bearer token) |

### External API

| Endpoint | Used By | Purpose |
|----------|---------|---------|
| `https://rest.uniprot.org/uniprotkb/:accession?format=json` | `UniprotProteinLabel` | Fetch protein metadata for hover tooltips |

---

## 10. Type Declarations

Three `.d.ts` files provide TypeScript type definitions for libraries that lack them:

- **`cytoscape-dagre.d.ts`** — Declares `cytoscape-dagre` as a Cytoscape extension
- **`cytoscape-svg.d.ts`** — Declares `cytoscape-svg` as a Cytoscape extension
- **`plotly.d.ts`** — Declares `plotly.js-dist-min` (re-exports from `plotly.js`) and `react-plotly.js/factory` (the `createPlotlyComponent` factory function with `PlotParams` interface)

---

## 11. Styling

### Framework

TailwindCSS 3.4 with **class-based dark mode** (`darkMode: 'class'` in `tailwind.config.js`).

### Dark Mode Architecture

The application supports light and dark themes with a user-toggleable switch:

1. **`ThemeContext.tsx`**: React context + provider that manages the current theme (`'light' | 'dark'`). Persists the user's preference to `localStorage` under the key `metagomics-theme`. On mount, reads the stored preference (defaults to `light`). The provider adds/removes the `dark` CSS class on `document.documentElement`, which activates all `dark:` Tailwind variants.

2. **`ThemeProvider`** wraps the entire app in `main.tsx` (above `BrowserRouter`), making `useTheme()` available to all components.

3. **`ThemeToggle`** (`components/ThemeToggle.tsx`): A pill-shaped toggle switch rendered in the Layout header. Uses Lucide `Sun` and `Moon` icons. The knob slides left/right with a CSS transition, and background icons indicate the inactive mode. Has an accessible `aria-label` that updates based on current state.

4. **Tailwind `dark:` variants**: Every UI element that has a light-mode color also specifies a `dark:` variant (e.g., `bg-white dark:bg-gray-900`, `text-gray-900 dark:text-gray-100`, `border-gray-200 dark:border-gray-700`).

5. **CSS custom properties** (`index.css`): Surface/text colors and DAG glow parameters are defined as CSS variables in `:root` and `.dark` selectors, consumed by `body` and global styles.

6. **Centralized JS color constants** (`utils/colors.ts`): All colors consumed programmatically by Cytoscape.js and Plotly.js are defined here (see Section 8.7). Components import light/dark color sets and select based on `useTheme()`.

### Color Scheme

| Element | Light Mode | Dark Mode |
|---------|-----------|-----------|
| **Page background** | `gray-50` (#f9fafb) | `gray-950` (#030712) |
| **Cards/panels** | `white` | `gray-900` |
| **Primary text** | `gray-900` | `gray-100` |
| **Secondary text** | `gray-600` | `gray-400` |
| **Borders** | `gray-200` | `gray-700` |
| **Primary accent** | `indigo-600` | `indigo-400` |
| **Links** | `indigo-600` | `indigo-400` |
| **Error** | `red-50` bg / `red-700` text | `red-900/30` bg / `red-400` text |
| **Status badges** | Semantic `100`-level bg / `800`-level text | Semantic `500/20` bg / `300`-level text with `ring-1` border |
| **Inputs/selects** | `white` bg / `gray-300` border | `gray-800` bg / `gray-600` border |
| **Tooltips** | `white` bg / `gray-300` border | `gray-800` bg / `gray-600` border with subtle indigo shadow |
| **Taxonomy palettes** | Medium-saturation colorblind-friendly (`DOMAIN_PALETTE_LIGHT`) | Vivid neon high-saturation (`DOMAIN_PALETTE_DARK`) |
| **Taxonomy text** | Auto (Plotly default) | Dark `#111827` for contrast on bright backgrounds |
| **GO DAG nodes** | Light tint backgrounds, gray borders | Indigo-bordered with glow shadow, vivid color ramp |
| **GO DAG edges** | `gray-300` | `gray-600`, highlighted ancestors in `indigo-300` |

### Layout

`max-w-7xl` centered container, responsive padding.

### Component Patterns

Cards use `bg-white dark:bg-gray-900 border border-gray-200 dark:border-gray-700 rounded-lg`. Consistent spacing via Tailwind utilities.

### Global Transitions

`index.css` applies `transition-property: background-color, border-color` to all elements for smooth dark mode transitions (150ms ease).

---

## 12. Build and Development

### Scripts (`package.json`)

| Script | Command | Purpose |
|--------|---------|---------|
| `dev` | `vite` | Start dev server with HMR and API proxy |
| `build` | `tsc && vite build` | Type-check then build for production |
| `lint` | `eslint . --ext ts,tsx` | Lint TypeScript/React files |
| `preview` | `vite preview` | Preview production build locally |
| `test` | `vitest run` | Run tests once |

### Vite Configuration (`vite.config.ts`)

- **Plugin**: `@vitejs/plugin-react`
- **Dev proxy**: `/api` → `http://localhost:8000` (to backend FastAPI server)
- **Test config**: jsdom environment, setup file `./src/test-setup.ts`, globals enabled

### TypeScript Configuration (`tsconfig.json`)

- **Target**: ES2020
- **Module**: ESNext with bundler resolution
- **Strict mode**: Enabled (`strict: true`)
- **Additional checks**: `noUnusedLocals`, `noUnusedParameters`, `noFallthroughCasesInSwitch`
- **JSX**: `react-jsx` (automatic runtime)
- **No emit**: `noEmit: true` (Vite handles bundling)

### Docker Build

The `Dockerfile` has a `frontend-builder` stage:
```dockerfile
FROM node:20-slim AS frontend-builder
WORKDIR /frontend
COPY frontend/package*.json ./
RUN npm ci
COPY frontend/ ./
RUN npm run build
```

The built output (`frontend/dist/`) is copied into the final Python image and served by FastAPI at `/assets/` (static) and `index.html` (SPA fallback).

---

## 13. Testing

### Philosophy

- **Utility functions have comprehensive tests**: Every parser in `utils/` has a dedicated test file covering happy paths, edge cases, empty input, malformed data, and realistic scenarios.
- **Component tests mock external dependencies**: `PeptideDetailsPane` tests mock `duckdb.ts` to avoid requiring a real WASM engine in tests.
- **Tests use `describe`/`it` blocks** with descriptive names following the `vitest` convention.
- **Testing Library for React components**: Tests use `render`, `screen`, `fireEvent`, `waitFor` from `@testing-library/react`.

### Running Tests

```bash
# In Docker (recommended — no local npm install needed):
docker build --target frontend-builder -t metagomics2-frontend-test .
docker run --rm metagomics2-frontend-test npx vitest run

# TypeScript type checking (also in Docker):
docker run --rm metagomics2-frontend-test npx tsc --noEmit
```

### Test Configuration

Defined in `vite.config.ts`:
```typescript
test: {
  environment: 'jsdom',
  setupFiles: ['./src/test-setup.ts'],
  globals: true,  // describe, it, expect available globally
}
```

Setup file (`test-setup.ts`): Imports `@testing-library/jest-dom` for DOM matchers like `toBeInTheDocument()`.

### Test Files

| Test File | Module Under Test | What It Verifies |
|-----------|-------------------|------------------|
| `utils/__tests__/csvParser.test.ts` | `csvParser.ts` | Simple fields, empty strings, quoted fields, escaped quotes, commas in quotes, newlines in quotes, realistic CSV lines |
| `utils/__tests__/taxonomyParser.test.ts` | `taxonomyParser.ts` | `parseTaxonomyCsv`: empty/header-only/single/multiple rows, quoted names, high precision decimals, short rows. `CANONICAL_RANKS`: membership. `filterCanonicalRanks`: canonical filtering, parent re-linking through non-canonical nodes, NCBI root normalization, realistic hierarchies. `filterByMaxRank`: all rank cutoffs. `validateCanonicalHierarchy`: valid hierarchies, rank gaps allowed, missing parents, same-rank parents, no-parent errors. `ensureStrictRankLayers`: no modification when consecutive, single/multiple placeholder insertion, placeholder naming, multiple children with same gap. |
| `utils/__tests__/goParser.test.ts` | `goParser.ts` | Empty/header-only/single/multiple rows, semicolon-delimited parents, empty parents, whitespace trimming, quoted names with commas, non-numeric defaults, short rows, high precision, Windows line endings |
| `utils/__tests__/comboParser.test.ts` | `comboParser.ts` | `parseComboCsv`: empty/header-only/single/multiple rows, empty parents, short rows. `comboRowsToTaxonNodes`: filtering by GO ID, field mapping, non-existent GO term. `comboRowsToGoTermNodes`: filtering by tax ID, field mapping, non-existent tax ID. |
| `components/__tests__/PeptideDetailsPane.test.tsx` | `PeptideDetailsPane` | Placeholder when no selection, loading state, correct peptide hierarchy rendering, "no peptides found" message, expand/collapse peptide rows, expand background protein to show annotated proteins, selection info in header (tax IDs and GO IDs). Mocks DuckDB via `vi.mock`. |
| `__tests__/ThemeContext.test.tsx` | `ThemeContext` | Default light theme, reads stored theme from localStorage, toggles light↔dark, persists to localStorage on toggle, ignores invalid stored values and defaults to light, adds/removes `dark` class on `document.documentElement` |
| `components/__tests__/ThemeToggle.test.tsx` | `ThemeToggle` | Renders accessible button with correct aria-label, toggles to dark mode on click, toggles back to light on second click, starts in dark mode when localStorage has dark preference |

### Writing New Tests

When adding new functionality:
1. **Utility functions**: Add tests in `utils/__tests__/test_<module>.test.ts`. Follow existing patterns — test empty input, single item, multiple items, edge cases, realistic data.
2. **Components**: Add tests in `components/__tests__/<Component>.test.tsx`. Mock external dependencies (API calls, DuckDB). Use Testing Library queries (`screen.getByText`, `screen.queryByText`).
3. **Run tests in Docker** before committing:
   ```bash
   docker build --target frontend-builder -t metagomics2-frontend-test .
   docker run --rm metagomics2-frontend-test npx vitest run
   docker run --rm metagomics2-frontend-test npx tsc --noEmit
   ```

### CI Pipeline (`.github/workflows/ci.yml`)

The `frontend-tests` job:
1. Builds the `frontend-builder` Docker target
2. Runs `npx tsc --noEmit` (TypeScript type checking)
3. Runs `npx vitest run` (unit tests)

---

## 14. Data Flow Summary

### Job Submission
```
NewJobPage → POST /api/jobs (FormData) → navigate to /job/:jobId
```

### Job Monitoring
```
JobPage → GET /api/jobs/:jobId (poll every 3s) → display status/progress → stop on completed/failed
```

### Taxonomy Visualization
```
GET taxonomy_nodes.csv → parseTaxonomyCsv → filterCanonicalRanks → validateCanonicalHierarchy
                          ↓
                     filterByMaxRank → ensureStrictRankLayers → TaxonomyChart (Plotly)
                          ↓ (if GO filter active)
GET go_taxonomy_combo.csv → parseComboCsv → comboRowsToTaxonNodes → same pipeline
```

### GO Visualization
```
GET go_terms.csv → parseGoTermsCsv → filter by namespace & minRatioTotal → GoDagViewer (Cytoscape)
                          ↓ (if taxon filter active)
GET go_taxonomy_combo.csv → parseComboCsv → comboRowsToGoTermNodes → same pipeline
```

### Peptide Drill-Down
```
User clicks chart node → selectedTaxIds / selectedGoId
    ↓
PeptideDetailsPane → registerMappingFile (DuckDB HTTP Parquet) → SQL query with WHERE clause
    ↓
Display: peptide → background_protein → annotated_protein (expandable tree)
```

---

## 15. Key Design Patterns

### Tooltip Positioning
Both `GoDagViewer` and `TaxonomyChart` use the same tooltip pattern:
1. On hover: capture rendered/client position, set tooltip state
2. `useLayoutEffect`: measure tooltip dimensions, flip horizontally/vertically if overflowing viewport
3. Tooltip is `position: fixed` with `pointer-events: none`

### Ref-Based Callback Stability
Both visualization components store `onNodeClick` in a `useRef` to avoid stale closures and unnecessary re-renders:
```typescript
const onNodeClickRef = useRef(onNodeClick)
onNodeClickRef.current = onNodeClick
```

### Memoized Plotly Props
`TaxonomyChart` heavily memoizes `plotlyData`, `plotlyLayout`, and `plotlyConfig` so that tooltip state changes don't trigger Plotly re-renders.

### DuckDB Singleton
`duckdb.ts` uses a module-level promise to ensure the WASM engine is initialized exactly once, even if `getDuckDB()` is called concurrently.

### UniProt Cache + In-Flight Deduplication
`uniprot.ts` maintains two module-level maps (`cache` and `inFlight`) to avoid duplicate API calls for the same accession.

### Theme-Aware Visualization Colors
Both `GoDagViewer` and `TaxonomyChart` consume the current theme via `useTheme()` and select the appropriate color set from the centralized `utils/colors.ts` module. Because Cytoscape and Plotly render to canvas/SVG (not the DOM), their colors cannot use CSS `dark:` variants and must be passed programmatically. The `isDark` boolean is included in `useMemo`/`useEffect` dependency arrays so visualizations re-render when the user toggles theme.

### Ancestor Path Highlighting (GO DAG)
On node hover, `GoDagViewer` performs a BFS traversal via `cy.incomers('edge')` to collect all ancestor nodes and edges up to root nodes, then applies the `highlighted` CSS class to the entire path. On mouseout, all highlights are cleared in a single `cy.nodes().removeClass()` / `cy.edges().removeClass()` call rather than tracking which elements were highlighted.

---

## 16. Backend Result Files Consumed by Frontend

| File | Format | Consumed By | How |
|------|--------|-------------|-----|
| `taxonomy_nodes.csv` | CSV | `TaxonomyPage`, `GoDagPage` (for taxon autocomplete) | `fetch` → `parseTaxonomyCsv` |
| `go_terms.csv` | CSV | `GoDagPage`, `TaxonomyPage` (for GO autocomplete) | `fetch` → `parseGoTermsCsv` |
| `go_taxonomy_combo.csv` | CSV | `GoDagPage` (taxon filter), `TaxonomyPage` (GO filter) | `fetch` → `parseComboCsv` → reshape |
| `coverage.csv` | CSV | Download only (JobPage) | Direct download link |
| `run_manifest.json` | JSON | Download only (JobPage) | Direct download link |
| `peptide_mapping.parquet` | Parquet | `PeptideDetailsPane` | DuckDB-WASM HTTP range requests |

---

## 17. Glossary

| Term | Definition |
|------|------------|
| **SPA** | Single-Page Application — the entire UI runs in one HTML page with client-side routing |
| **Canonical ranks** | root, domain, kingdom, phylum, class, order, family, genus, species — the 9 standard taxonomic levels displayed in charts |
| **Placeholder node** | A synthetic node inserted by `ensureStrictRankLayers` to fill rank gaps in the taxonomy hierarchy, ensuring correct chart ring depth |
| **Namespace** | GO ontology division: biological_process, cellular_component, or molecular_function |
| **DuckDB-WASM** | WebAssembly build of DuckDB that runs SQL queries in the browser without a server |
| **Combo data** | The `go_taxonomy_combo.csv` cross-tabulation of taxonomy × GO with quantities and fractions for each combination |
| **fractionOfTaxon** | For a (taxon, GO) combo: `quantity / taxon_total_quantity` — what fraction of this taxon's quantity is attributed to this GO term |
| **fractionOfGo** | For a (taxon, GO) combo: `quantity / go_total_quantity` — what fraction of this GO term's quantity is attributed to this taxon |
