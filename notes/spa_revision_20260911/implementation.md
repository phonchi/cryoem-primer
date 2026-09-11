# SPA revision — 2026-09-11

Audience: Traditional Chinese; undergraduates with basic statistics, calculus and linear algebra.
Scope: Expand SPA only; keep general-image chapters 01–04 and synthetic chapter 07, their py/ipynb/text/outputs/assets byte-identical. Protected hashes: protected_files.json. Do not globally sync notebooks or change their tests to remove existing comprehension checks.
SPA order: 05_background, 05_image_formation, 05_statistical_inference, 06_workflow, 06_alignment_classification, 06_dimension_reduction, 06_reconstruction_validation, 06_heterogeneity, 06_resolution_validation.
Root owns shared toc/bib/integration; preserve old filenames and anchors. Sources remain local; publish only bibliographic metadata and original teaching prose. No push/deploy authorized yet.
Teaching: motivation, intuition, core equations, worked examples in main text; long derivations in dropdowns. No comprehension question sections in revised pages. Each chapter has annotated full reading list, links and exact reading locators. Existing five protected chapters retain their checks.
Notation: image-plane-to-volume R; P_R V(u)=integral V(R[u,z]) dz; Fourier slice is Vhat(R[kx,ky,0]). y_i=A_i(z_i)V+noise. Explain likelihood/prior assumptions; distinguish denoising, classification and reconstruction.
Collaboration: Claude writes background, formation, workflow, alignment, reconstruction, heterogeneity, validation. Codex writes inference, dimension reduction, navigation, appendix, figures and integration. Native literature agent owns bibliography, document cache, reference routing and source digests.
Validation: protected-file SHA checks, scientific numerical examples, site/source/reference tests, cached whole-book build and browser inspection. Save outputs/logs here or logs/spa_revision_20260911. Do not introduce paid services or alter existing dependencies.

User steering: replace student-facing 姿態 with 角度; use 角度與位移 for joint rotation/translation. Apply speak-human-tw proofreading without changing protected chapters or scientific content.

Final user authorization: after scientific and language corrections, dispatch multiple independent readers to read the entire book and assess continuity; fix the findings within the allowed SPA/navigation scope, then commit and push and verify GitHub Pages deployment. This supersedes the earlier publication-approval-pending statement. Protected five chapters/assets remain unchanged.
