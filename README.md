# Costim Counter-Screen — Built with Claude (Life Sciences · Researcher Track)

Max Goulazian — nomination: **4-1BB (TNFRSF9) + CD27**

This folder is the **complete, self-contained web root**. Static files only — no server,
no build step. Open `index.html` (it redirects to `Costim - Start Here.dc.html`).

## Push to GitHub Pages

You're already inside the web root, so just:

```bash
git init
git add -A
git commit -m "Costim counter-screen — Built with Claude"
git branch -M main
git remote add origin https://github.com/<you>/<repo>.git
git push -u origin main
```

Then **Settings → Pages → Deploy from a branch → `main` / `root` → Save**.
Live at `https://<you>.github.io/<repo>/`.

Every file here pushes as normal git (largest ≈ 3–4 MB, far under the 100 MB limit) —
**no Git LFS needed**. `.nojekyll` is included so Pages serves the `.dc.html` pages and
asset folders verbatim.

> Use the **git command line**, not GitHub's drag-and-drop uploader — the web uploader
> caps at ~25 MB / 100 files and can't take the full asset set.

## The demo video

The landing page has an open **16:9 video slot** pointing at **`demo.mp4`**. Record the
walkthrough, name it `demo.mp4`, and drop it in this folder next to `index.html` — it fills
the slot automatically. Until then the slot shows the binder still as a poster.

## What's inside
- **Landing** — `Costim - Start Here.dc.html` (+ `index.html` entry)
- **Deck** — `Costim Counter-Screen Deck.dc.html`
- **Interactive tools** — Landscape, PBPK Runner (spatial atlas: 23 antibodies × organ),
  QSP Model Explorer, GRN Networks, Binder Structure (3D), Binder Diffusion
- **Documents** — Figures, Writeup, QSP Model Documentation, Built with Claude
- **Notebooks** — reproduction `.ipynb` (rendered viewers + downloadable files)
- **Data** — `*_data.js`; **images** — `spatial_resolved/`, `spatial2/`, `grn/`,
  `deck_assets/`, `timelapse/`

Needs internet at runtime: Google Fonts + the 3Dmol structure viewer load from CDNs.
