# Push to GitHub + go live on Pages — step by step

A fully static site (no build step). Goal: same behavior as the preview, served at
`https://<you>.github.io/<repo>/`.

## 0. One-time prerequisites
- A GitHub account.
- Git installed locally (`git --version` to check).
- This folder unzipped somewhere you can `cd` into. **This folder IS the web root** — the
  files inside (`index.html`, the `.dc.html` pages, `sweep/`, etc.) go at the repo root,
  not inside a subfolder.

## 1. Create an empty repo on GitHub
1. GitHub → top-right **+** → **New repository**.
2. Name it (e.g. `costim-counterscreen`).
3. **Public** (required for free GitHub Pages).
4. Do **NOT** add a README, .gitignore, or license (keep it empty — avoids a merge conflict).
5. **Create repository.** Leave the page open; you'll need the URL.

## 2. Push this folder
In a terminal, from inside the unzipped `github_site/` folder:

```bash
cd path/to/github_site        # the folder with index.html in it
git init
git add -A
git commit -m "Costim counter-screen — Built with Claude"
git branch -M main
git remote add origin https://github.com/<you>/<repo>.git
git push -u origin main
```

- No Git LFS needed — every file is well under 100 MB.
- If it asks for a password, use a **Personal Access Token** (GitHub no longer accepts your
  account password on the CLI): GitHub → Settings → Developer settings → Tokens. Use a
  fine-grained or classic token with `repo` scope, paste it as the password.
- Filenames with spaces (e.g. `Costim - Start Here.dc.html`) are fine — git and Pages handle them.

## 3. Turn on GitHub Pages
1. Your repo → **Settings** → **Pages** (left sidebar).
2. **Build and deployment → Source: Deploy from a branch.**
3. **Branch: `main`**, **Folder: `/ (root)`** → **Save.**
4. Wait ~1–2 minutes. The page will show your live URL: `https://<you>.github.io/<repo>/`.

`index.html` auto-redirects to `Costim - Start Here.dc.html`, so the landing page loads first
and every deck / figure / tool / notebook links from there.

## 4. Confirm it's fully functional
Open the live URL and check:
- Landing loads, then click into the **Deck**, **Figures**, **PBPK Runner**, **QSP Explorer**,
  **Binder Structure (3D)**, **Diffusion**, **Writeup**, notebooks.
- The 3D binder viewer and Google fonts load from CDNs → needs internet (they will not work offline).
- Spatial atlas + sweep images load on scroll (lazy-loaded — first view of each may take a moment).

`.nojekyll` is already included so Pages serves the `.dc.html` files and asset folders verbatim
(without it, GitHub's Jekyll step can drop files).

## 5. Add your demo video (optional)
The landing page has a 16:9 slot wired to `demo.mp4`.
```bash
cp /path/to/your-walkthrough.mp4 demo.mp4    # into the repo root, next to index.html
git add demo.mp4 && git commit -m "Add demo video" && git push
```
It fills the slot automatically. (If the mp4 is >100 MB, either compress it or track it with
Git LFS — but LFS files do NOT play on GitHub Pages, so prefer keeping it under 100 MB.)

## 6. Updating later
Any change: `git add -A && git commit -m "..." && git push`. Pages redeploys in ~1 min.

## Serve as the landing page of an EXISTING submission repo

Recommended: a dedicated `gh-pages` branch. The site lives at that branch's root (so all
relative links work unchanged), and your submission files + LFS on `main` are untouched.

```bash
git clone https://github.com/<you>/<repo>.git
cd <repo>
git checkout --orphan gh-pages     # new empty branch; main untouched
git rm -rf .                        # clears working tree on gh-pages only
cp -R /path/to/github_site/. .      # copy the unzipped bundle contents to root
git add -A
git commit -m "Site landing page"
git push -u origin gh-pages
```
Then Settings → Pages → Deploy from a branch → Branch **`gh-pages`** / **`/ (root)`** → Save.
Live at `https://<you>.github.io/<repo>/`; submission files stay on `main`.

**LFS caveat:** if your submission repo's `.gitattributes` LFS-tracks `*.png` globally, that rule
would catch the site's spatial/sweep PNGs — and GitHub Pages does NOT serve LFS bytes, so images
break. The `gh-pages` orphan branch starts with no LFS rules, so the site's PNGs commit as normal
git and serve fine. Alternative `/docs` route (Pages → main / `/docs`): first add
`docs/** -filter -diff -merge` to `.gitattributes` so the site images aren't LFS-tracked.

## Troubleshooting
- **404 at the URL:** Pages source must be `main` / **root** (not `/docs`); wait a minute after Save.
- **Blank pages / missing files:** confirm `.nojekyll` is in the repo root (it is in this bundle).
- **Images missing:** don't add site assets to `.gitattributes`/LFS — Pages won't serve LFS bytes.
- **"src refspec main does not match":** you committed before `git branch -M main`; re-run step 2 in order.
