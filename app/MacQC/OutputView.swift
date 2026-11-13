//
//  OutputView.swift
//  MacQC
//
//  Output tab: the output file of the job, streamed while it runs (H₂, LiH, H₂O), or the
//  FORTRAN-style report of the Szabo–Ostlund calculation (HeH⁺). A read-only NSTextView keeps
//  long files fast; new text is appended, the view follows the end while "Follow output" is on,
//  the text is selectable, and the search field finds the next or previous occurrence.
//

import SwiftUI
import AppKit

struct OutputView: View {
    @Bindable var model: AppModel
    @State private var controller = OutputTextController()
    @State private var query = ""
    @State private var follow = true
    @State private var searchMessage = ""

    /// The text shown and when it was replaced rather than extended.
    private var content: (text: String, generation: Int) {
        if model.isGeneral || model.result?.isGeneral == true {
            if !model.outputText.isEmpty { return (model.outputText, model.outputGeneration) }
            if let o = model.result?.output { return (o, o.hashValue) }
            return ("", 0)
        }
        let report = model.result?.report ?? ""
        return (report, report.hashValue)
    }

    var body: some View {
        let c = content
        VStack(spacing: 0) {
            toolbar(c.text)
            Divider()
            if c.text.isEmpty {
                VStack(spacing: 14) {
                    if let message = model.errorMessage {
                        ErrorBanner(message: message)
                            .frame(maxWidth: 560)
                    }
                    ContentUnavailableView(
                        model.isRunning ? "Waiting for output…" : "No output yet",
                        systemImage: "doc.plaintext",
                        description: Text("The output file of the job appears here while it runs: input, molecule, basis set, every SCF, coupled-cluster and optimisation iteration, frequencies, thermochemistry and timings. File ▸ Export Results as Text saves the same file.")
                    )
                }
                .padding()
                .frame(maxWidth: .infinity, maxHeight: .infinity)
            } else {
                OutputTextView(text: c.text, generation: c.generation, follow: follow, controller: controller)
            }
        }
    }

    private func toolbar(_ text: String) -> some View {
        HStack(spacing: 10) {
            Image(systemName: "magnifyingglass")
                .foregroundStyle(.secondary)
            TextField("Search the output", text: $query)
                .textFieldStyle(.roundedBorder)
                .frame(maxWidth: 260)
                .onSubmit { find(backwards: false) }
                .accessibilityLabel("Search the output")
            Button {
                find(backwards: true)
            } label: {
                Image(systemName: "chevron.up")
            }
            .help("Previous match")
            .accessibilityLabel("Previous match")
            .disabled(query.isEmpty || text.isEmpty)
            Button {
                find(backwards: false)
            } label: {
                Image(systemName: "chevron.down")
            }
            .help("Next match (Return)")
            .accessibilityLabel("Next match")
            .disabled(query.isEmpty || text.isEmpty)
            if !searchMessage.isEmpty {
                Text(searchMessage)
                    .font(.caption)
                    .foregroundStyle(.secondary)
            }
            Spacer()
            if model.isRunning {
                ProgressView()
                    .controlSize(.small)
                Text("Running")
                    .font(.caption)
                    .foregroundStyle(.secondary)
            }
            Toggle("Follow output", isOn: $follow)
                .toggleStyle(.checkbox)
                .help("Keep the end of the output in view while the job runs.")
            Text(lineCount(text))
                .font(.caption.monospacedDigit())
                .foregroundStyle(.secondary)
        }
        .padding(.horizontal, 14)
        .padding(.vertical, 8)
    }

    private func lineCount(_ text: String) -> String {
        guard !text.isEmpty else { return "" }
        var n = 0
        for b in text.utf8 where b == 10 { n += 1 }
        return "\(n) lines"
    }

    private func find(backwards: Bool) {
        guard !query.isEmpty else { return }
        follow = false   // searching stops following the end
        let result = controller.find(query, backwards: backwards)
        switch result {
        case .found(let index, let total): searchMessage = "\(index) of \(total)"
        case .notFound: searchMessage = "not found"
        }
    }
}

/// Holds the text view so the search field can drive it.
@MainActor
final class OutputTextController {
    weak var textView: NSTextView?

    enum FindResult {
        case found(Int, Int)
        case notFound
    }

    /// Selects and shows the next (or previous) case-insensitive occurrence, wrapping around.
    func find(_ query: String, backwards: Bool) -> FindResult {
        guard let tv = textView, let storage = tv.textStorage, !query.isEmpty else { return .notFound }
        let text = storage.string as NSString
        let length = text.length
        let selected = tv.selectedRange()
        var range: NSRange
        if backwards {
            let end = selected.location
            range = text.range(of: query, options: [.caseInsensitive, .backwards], range: NSRange(location: 0, length: end))
            if range.location == NSNotFound {
                range = text.range(of: query, options: [.caseInsensitive, .backwards], range: NSRange(location: 0, length: length))
            }
        } else {
            let start = min(selected.location + selected.length, length)
            range = text.range(of: query, options: .caseInsensitive, range: NSRange(location: start, length: length - start))
            if range.location == NSNotFound {
                range = text.range(of: query, options: .caseInsensitive, range: NSRange(location: 0, length: length))
            }
        }
        guard range.location != NSNotFound else { return .notFound }
        tv.setSelectedRange(range)
        tv.scrollRangeToVisible(range)
        tv.showFindIndicator(for: range)
        // position of this match among all matches
        var count = 0, index = 0
        var searchRange = NSRange(location: 0, length: length)
        while true {
            let r = text.range(of: query, options: .caseInsensitive, range: searchRange)
            if r.location == NSNotFound { break }
            count += 1
            if r.location == range.location { index = count }
            let next = r.location + max(r.length, 1)
            if next >= length { break }
            searchRange = NSRange(location: next, length: length - next)
        }
        return .found(index, count)
    }
}

/// Read-only, selectable, monospaced text that grows by appending.
struct OutputTextView: NSViewRepresentable {
    var text: String
    var generation: Int
    var follow: Bool
    var controller: OutputTextController

    final class Coordinator {
        var generation = Int.min
        var length = 0   // UTF-16 length shown
    }

    func makeCoordinator() -> Coordinator { Coordinator() }

    func makeNSView(context: Context) -> NSScrollView {
        let scroll = NSTextView.scrollableTextView()
        scroll.hasHorizontalScroller = true
        scroll.autohidesScrollers = true
        if let tv = scroll.documentView as? NSTextView {
            tv.isEditable = false
            tv.isSelectable = true
            tv.isRichText = false
            tv.usesFindBar = true
            tv.isIncrementalSearchingEnabled = true
            tv.font = NSFont.monospacedSystemFont(ofSize: 11.5, weight: .regular)
            tv.textColor = .textColor
            tv.backgroundColor = .textBackgroundColor
            tv.textContainerInset = NSSize(width: 10, height: 10)
            // no line wrapping: the output is a fixed-width table layout
            tv.isHorizontallyResizable = true
            tv.maxSize = NSSize(width: CGFloat.greatestFiniteMagnitude, height: CGFloat.greatestFiniteMagnitude)
            tv.textContainer?.widthTracksTextView = false
            tv.textContainer?.containerSize = NSSize(width: CGFloat.greatestFiniteMagnitude,
                                                     height: CGFloat.greatestFiniteMagnitude)
            tv.setAccessibilityLabel("Output file")
            controller.textView = tv
        }
        return scroll
    }

    func updateNSView(_ scroll: NSScrollView, context: Context) {
        guard let tv = scroll.documentView as? NSTextView, let storage = tv.textStorage else { return }
        controller.textView = tv
        let c = context.coordinator
        let ns = text as NSString
        let attributes: [NSAttributedString.Key: Any] = [
            .font: NSFont.monospacedSystemFont(ofSize: 11.5, weight: .regular),
            .foregroundColor: NSColor.textColor,
        ]
        if c.generation != generation || ns.length < c.length {
            storage.setAttributedString(NSAttributedString(string: text, attributes: attributes))
            c.generation = generation
            c.length = ns.length
            if follow { tv.scrollToEndOfDocument(nil) }
        } else if ns.length > c.length {
            let piece = ns.substring(from: c.length)
            storage.append(NSAttributedString(string: piece, attributes: attributes))
            c.length = ns.length
            if follow { tv.scrollToEndOfDocument(nil) }
        }
    }
}
