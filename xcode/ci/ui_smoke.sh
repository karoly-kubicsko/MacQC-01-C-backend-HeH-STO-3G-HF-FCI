#!/bin/bash
# UI smoke test for the built MacQC.app on a GitHub macOS runner.
# Launches the app, screenshots each tab, dumps the visible text of every tab to the log,
# runs the potential-curve scan and checks the Edit/File menu commands.
# Never fails the job: it only gathers evidence (screenshots + text) for review.
#
# Usage: ui_smoke.sh path/to/MacQC.app output-dir

APP="$1"
OUT="$2"
mkdir -p "$OUT"

# Reads the text of the front window; optionally clicks a tab/button first.
cat > "$OUT/tab.applescript" <<'APPLESCRIPT'
on run argv
    set targetName to item 1 of argv
    tell application "System Events"
        tell process "MacQC"
            set frontmost to true
            if targetName is not "" then
                set els to entire contents of window 1
                repeat with e in els
                    try
                        if (name of e) is targetName then
                            if (role of e) is in {"AXRadioButton", "AXButton", "AXTab"} then
                                click e
                                exit repeat
                            end if
                        end if
                    end try
                end repeat
            end if
            delay 2
            set resultText to ""
            set els to entire contents of window 1
            repeat with e in els
                try
                    if (role of e) is "AXStaticText" then
                        set resultText to resultText & (value of e as text) & linefeed
                    end if
                end try
            end repeat
            return resultText
        end tell
    end tell
end run
APPLESCRIPT

shot() { screencapture -x "$OUT/$1.png" || echo "screencapture failed for $1"; }

dump_tab() {   # $1 = file prefix, $2 = element to click ("" = none)
    echo "::group::UI text: $1"
    timeout 120 osascript "$OUT/tab.applescript" "$2" | tee "$OUT/$1.txt"
    echo "::endgroup::"
    shot "$1"
}

echo "== launching $APP"
open -n "$APP"
sleep 15
pgrep -x MacQC >/dev/null && echo "MacQC process is running" || echo "MacQC process NOT running"

shot "00-launch"
dump_tab "01-summary" ""
dump_tab "02-hartree-fock" "Hartree–Fock"
dump_tab "03-full-ci" "Full CI"
dump_tab "04-potential-curve" "Potential curve"

echo "== running the potential-curve scan"
timeout 120 osascript "$OUT/tab.applescript" "Run scan" > /dev/null
sleep 10
dump_tab "05-curve-after-scan" ""

echo "== menu commands"
echo "::group::menu items"
timeout 60 osascript <<'APPLESCRIPT' | tee "$OUT/menus.txt"
tell application "System Events"
    tell process "MacQC"
        set out to ""
        repeat with m in {"File", "Edit"}
            try
                set names to name of every menu item of menu 1 of menu bar item (m as text) of menu bar 1
                set out to out & m & ": " & (names as text) & linefeed
            end try
        end repeat
        return out
    end tell
end tell
APPLESCRIPT
echo "::endgroup::"

echo "== Edit > Copy Summary as Text"
pbcopy < /dev/null
timeout 60 osascript <<'APPLESCRIPT' || echo "could not click the menu item"
tell application "System Events"
    tell process "MacQC"
        click menu item "Copy Summary as Text" of menu 1 of menu bar item "Edit" of menu bar 1
    end tell
end tell
APPLESCRIPT
sleep 1
echo "::group::clipboard after Copy Summary as Text"
pbpaste | tee "$OUT/clipboard-summary.txt"
echo "::endgroup::"

pkill -x MacQC || true
exit 0
