# GitHub and Zenodo

Use `mhv_polya` as the repository name and this folder as the repository root. The initial GitHub upload can be made now; complete the publication metadata before creating the archived release.

## 1. Upload the repository

Create an empty repository named `mhv_polya` on GitHub. Leave the automatic README, .gitignore and license options unselected, because the first commit supplies the local files. Copy the repository's HTTPS URL. These steps follow the [GitHub import guide](https://docs.github.com/en/migrations/importing-source-code/using-the-command-line-to-import-source-code/adding-locally-hosted-code-to-github).

Run the following commands in Terminal. Replace `OWNER` with your GitHub username or organization:

```bash
cd /path/to/mhv_polya
git init -b main
git add .
git diff --cached --stat
git commit -m "Add MHV poly(A) figure reproduction"
git remote add origin https://github.com/OWNER/mhv_polya.git
git push -u origin main
```

The commit includes scripts, README/docs, the figure preview, method settings, table schemas and expected TSV summaries. `.gitignore` excludes sequencing data, generated runs, Parquet/NPY intermediates, local archives, caches and local credentials. GitHub authentication must be configured for the push.

For subsequent edits:

```bash
git add .
git diff --cached --stat
git commit -m "Update analysis documentation"
git push
```

## 2. Complete the publication version

The release needs the following information:

| Item | What to supply | Where it goes |
|---|---|---|
| Software citation | Software authors in order, version and GitHub repository URL | Root `CITATION.cff`, using [the template](../config/CITATION.cff.template) |
| Software license | The authors' selected license, with the applicable copyright details | Root `LICENSE`; matching license identifier in `CITATION.cff` |
| Input data | Accession or DOI for the 36 required input files | README's **Input data accession** entry |
| Archived code | DOI for the released code version, assigned after Zenodo processes the release | README's **Code DOI** entry, `CITATION.cff` and manuscript |

Fill in the citation template, then save it as `CITATION.cff` in the repository root. GitHub uses this file for its [Cite this repository](https://docs.github.com/en/repositories/managing-your-repositorys-settings-and-features/customizing-your-repository/about-citation-files) feature. Add the selected `LICENSE`. The code DOI can be added after the first release is archived.

Use an existing data accession if it provides all 36 required inputs; otherwise deposit the remaining inputs and link their record. Include the five viral BAM indices if releasing the advanced workflow's input set. The code DOI identifies the software archive; the data accession identifies the study inputs. Record the data accession in the README and GitHub release description.

The [validation notes](VALIDATION.md) describe the remaining installation check.

In the analysis environment, run:

```bash
python scripts/check_bundle.py --release
python scripts/build_archive.py --output /path/to/mhv_polya_release.zip
```

Commit and push the completed metadata before creating the release.

## 3. Obtain the code DOI

Connect GitHub to Zenodo and [enable `mhv_polya`](https://help.zenodo.org/docs/github/enable-repository/). After enabling the repository, create a GitHub release from the manuscript's code version, for example `v1.0.0`. Zenodo processes that release and provides its DOI. Follow the [Zenodo release guide](https://help.zenodo.org/docs/github/archive-software/github-upload/).

Check the Zenodo record's title, authors, version and license. Add the version DOI to the manuscript and citation metadata. Release notes should identify the analysis version, data accession and validation report. GitHub displays the root [CITATION.cff](https://docs.github.com/en/repositories/managing-your-repositorys-settings-and-features/customizing-your-repository/about-citation-files) as the repository citation.

Replace the two pending entries in the README with the actual code DOI and input-data accession. Commit those documentation updates to the default branch. The already archived release remains the snapshot of its tagged commit.
