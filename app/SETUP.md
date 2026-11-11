# ΨLab app – Xcode setup

The SwiftUI sources live in `app/MacQC/`. No `.xcodeproj` is committed; you create it once with
the steps below (Xcode 26, macOS 14 or later as deployment target). The app is a universal binary:
it runs natively on Apple Silicon and on Intel Macs.

Paths below assume the repository root is `<repo>` and the Xcode project is created in
`<repo>/xcode/`, so that `$(SRCROOT)` = `<repo>/xcode`.

## Shortcut: generate the project instead of creating it by hand

`xcode/project.yml` is an [XcodeGen](https://github.com/yonaskolb/XcodeGen) spec that encodes every
setting below (sources, bridging header, header search paths, C++17, macOS 14.0, `-ffp-contract=off`,
sandbox entitlement, universal build):

```bash
brew install xcodegen
cd xcode && xcodegen generate     # creates xcode/PsiLab.xcodeproj; open it and Run
```

Then skip to step 5. The generated project is git-ignored. GitHub Actions
(`.github/workflows/app.yml`) uses the same spec to compile-check the app on a macOS runner and to
publish a universal, ad-hoc signed `ΨLab.app` as the workflow artifact `PsiLab-app`.
To run it: download `PsiLab-app` from the run's *Artifacts* section; GitHub delivers it as a zip,
which contains `ΨLab.app.zip`; unzip that too to get `ΨLab.app`. It is not notarized, so either
right-click ▸ Open once, or run `xattr -dr com.apple.quarantine ΨLab.app` (type the command up to
the space, then drag the app onto the Terminal window to insert its path).

## 1. Create the project (manual route)

1. Xcode ▸ *File ▸ New ▸ Project…* ▸ **macOS** ▸ **App** ▸ Next.
2. Product Name **PsiLab**, Interface **SwiftUI**, Language **Swift**, Testing System *None*,
   Storage *None*. Next.
3. Choose (create) the folder `<repo>/xcode` and **uncheck** "Create Git repository". Create.
4. In the Project navigator delete the template's `PsiLabApp.swift` and `ContentView.swift`
   (*Move to Trash*). Keep `Assets.xcassets` (and the template's entitlements file).

## 2. Add the sources

Use *File ▸ Add Files to "PsiLab"…* (or drag from Finder). Always: **"Copy items if needed" off**,
"Create groups", and the PsiLab target ticked.

| What | Files | Target membership |
| --- | --- | --- |
| App | everything in `app/MacQC/` (12 `.swift` files and `MacQC-Bridging-Header.h`) | `.swift`: PsiLab; header: none needed |
| Backend sources | `backend/src/linalg.cpp`, `integrals.cpp`, `report.cpp`, `scf.cpp`, `fci.cpp`, `calculation.cpp`, `api.cpp` | PsiLab |
| Backend headers | `backend/include/macqc/*.h` | (headers only) |
| JSON library | `third_party/nlohmann/json.hpp` | (header only) |

**Do not add** `backend/cli/main.cpp`, anything in `backend/tests/`, or `reference/` (the CLI and
the test program have their own `main`, which would clash with the app's).

Check *PsiLab target ▸ Build Phases ▸ Compile Sources*: exactly 12 Swift files + 7 C++ files.

## 3. Build settings

Select the project ▸ target **PsiLab** ▸ *Build Settings* ▸ "All" and "Combined", and set:

| Setting | Value |
| --- | --- |
| macOS Deployment Target (`MACOSX_DEPLOYMENT_TARGET`) | **14.0** |
| Product Name (`PRODUCT_NAME`) | `ΨLab` (the visible app name) |
| Product Module Name (`PRODUCT_MODULE_NAME`) | `PsiLab` (keeps the Swift module name ASCII) |
| Objective-C Bridging Header (`SWIFT_OBJC_BRIDGING_HEADER`) | `$(SRCROOT)/../app/MacQC/MacQC-Bridging-Header.h` |
| Header Search Paths (`HEADER_SEARCH_PATHS`) | `$(SRCROOT)/../backend/include` and `$(SRCROOT)/../third_party` |
| C++ Language Dialect (`CLANG_CXX_LANGUAGE_STANDARD`) | **C++17** (`c++17`; GNU++17 also works) |
| C++ Standard Library (`CLANG_CXX_LIBRARY`) | libc++ (default) |
| Other C++ Flags (`OTHER_CPLUSPLUSFLAGS`) | `-ffp-contract=off` (keeps results bit-identical to the reference; see backend README) |
| Architectures (`ARCHS`) | Standard Architectures (Apple Silicon, Intel) – the default |
| Build Active Architecture Only, Release | No (for a universal release build) |

Leave *C++ and Objective-C Interoperability* at "C / Objective-C" (the bridging header only
includes the plain C API). The Swift language mode may be Swift 6 or Swift 5; the code is written
to compile in both, and with "Default Actor Isolation" set either to `MainActor` (Xcode 26 template
default) or `nonisolated`.

## 4. Signing & Capabilities

- *Signing*: "Sign to Run Locally" or your personal team is enough for "My Mac".
- *App Sandbox* (enabled by the template): under **File Access ▸ User Selected File** choose
  **Read/Write** – without it, File ▸ Export Results as JSON… cannot save. No network or other
  entitlements are needed.

## 5. Build and run

1. Choose the scheme **PsiLab** and the destination **My Mac** in the toolbar.
2. *Product ▸ Run* (⌘R in Xcode). The app computes the Szabo–Ostlund default case on launch.
3. Expected for the defaults (STO-3G, R = 1.4632 bohr, ζ = 2.0925 / 1.24):
   HF total −2.86066199 E_h (6 iterations), FCI total −2.88071185 E_h, correlation −0.02004986 E_h.

Inside the app, ⌘R runs the calculation, File ▸ *Export Results as JSON…* (⇧⌘E) saves the result
and Edit ▸ *Copy Summary as Text* (⇧⌘C) copies a plain-text summary.

## Liquid Glass

The glass effects (`glassEffect`, `GlassEffectContainer`, `.glass`, `.glassProminent`) are used only
when running on macOS 26 (`if #available(macOS 26, *)`); on macOS 14–15 the app falls back to
`.regularMaterial` and the standard bordered button styles. Compiling requires the **macOS 26 SDK
(Xcode 26)**, even for a deployment target of 14.0.

## Troubleshooting

| Symptom | Fix |
| --- | --- |
| `'macqc/macqc_api.h' file not found` | Header Search Paths missing/wrong (`$(SRCROOT)` is the folder with `PsiLab.xcodeproj`). |
| `Cannot find 'macqc_run' in scope` | Bridging header path not set, or the path does not point to `app/MacQC/MacQC-Bridging-Header.h`. |
| `'nlohmann/json.hpp' file not found` | Add `$(SRCROOT)/../third_party` to Header Search Paths. |
| `duplicate symbol '_main'` | `backend/cli/main.cpp` or the test program was added to the target; remove it. |
| `Cannot find 'glassEffect' / 'glassProminent'` | Xcode/SDK older than 26. |
| Export does nothing / error on saving | App Sandbox ▸ User Selected File is not Read/Write. |
| "@main" or "invalid redeclaration" errors | The template's `PsiLabApp.swift` / `ContentView.swift` were not deleted. |
