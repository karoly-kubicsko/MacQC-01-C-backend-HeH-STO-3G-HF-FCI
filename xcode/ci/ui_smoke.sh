#!/bin/bash
# UI smoke test for the built ΨLab.app on a GitHub macOS runner.
# Launches the app, screenshots each tab, dumps the visible text of every tab to the log,
# runs the potential-curve scan and checks the Edit/File menu commands; then relaunches with
# launch arguments for H2O (RHF), LiH (Full CI), an H2 potential curve, the refused Full CI
# of H2O, a long H2O calculation (warning) and the coupled-cluster methods of ΨLab 4.0
# (open method menu, H2O CCSD(T), LiH CCSDT, refused H2O CCSDTQ).
# Never fails the job: it only gathers evidence (screenshots + text) for review.
#
# Usage: ui_smoke.sh path/to/ΨLab.app output-dir

APP="$1"
OUT="$2"
mkdir -p "$OUT"

# Reads the text of the front window; optionally clicks a tab/button first.
cat > "$OUT/tab.applescript" <<'APPLESCRIPT'
on run argv
    set targetName to item 1 of argv
    set logText to ""
    tell application "System Events"
        tell process "ΨLab"
            set frontmost to true
            if targetName is "__TREE__" then
                -- accessibility tree (everything except static text), to find how tabs are exposed
                set els to entire contents of window 1
                repeat with e in els
                    try
                        set r to role of e
                        if r is not "AXStaticText" then
                            set n to ""
                            set d to ""
                            try
                                set n to name of e
                            end try
                            try
                                set d to description of e
                            end try
                            set logText to logText & r & " | " & n & " | " & d & linefeed
                        end if
                    end try
                end repeat
                return logText
            end if
            if targetName is "__SCAN__" then
                -- the Run scan button has no accessible name: click the first AXButton
                -- that follows the "Step (Å)" pop-up in the tree
                set els to entire contents of window 1
                set seenStep to false
                repeat with e in els
                    try
                        set r to role of e
                        set n to ""
                        try
                            set n to name of e
                        end try
                        if r is "AXPopUpButton" and n is "Step (Å)" then
                            set seenStep to true
                        else if seenStep and r is "AXButton" then
                            click e
                            return "CLICKED button after Step pop-up"
                        end if
                    end try
                end repeat
                return "NOT FOUND: button after Step pop-up"
            end if
            if targetName is not "" then
                set clicked to false
                set els to entire contents of window 1
                repeat with e in els
                    if clicked then exit repeat
                    try
                        set r to role of e
                        if r is not "AXStaticText" and r is not "AXImage" and r is not "AXGroup" and r is not "AXScrollArea" then
                            set n to ""
                            set d to ""
                            set t to ""
                            try
                                set n to name of e
                            end try
                            try
                                set d to description of e
                            end try
                            try
                                set t to title of e
                            end try
                            if n is targetName or d is targetName or t is targetName then
                                click e
                                set clicked to true
                                set logText to logText & "CLICKED " & r & " | " & n & " | " & d & linefeed
                            end if
                        end if
                    end try
                end repeat
                if not clicked then set logText to logText & "NOT FOUND: " & targetName & linefeed
            end if
            delay 2
            set els to entire contents of window 1
            repeat with e in els
                try
                    if (role of e) is "AXStaticText" then
                        set logText to logText & (value of e as text) & linefeed
                    end if
                end try
            end repeat
            return logText
        end tell
    end tell
end run
APPLESCRIPT

# macOS has no `timeout`: run a command and kill it after N seconds.
with_timeout() {
    local secs="$1"; shift
    "$@" &
    local pid=$!
    ( sleep "$secs"; kill "$pid" 2>/dev/null ) > /dev/null 2>&1 &
    local watchdog=$!
    wait "$pid"
    local rc=$?
    kill "$watchdog" 2>/dev/null
    return $rc
}

# Screenshot, and print it at once into the log as a downscaled base64 JPEG, so the images can
# be read through the GitHub API without downloading the artifact (and survive a timeout).
# Decode: the line between "IMG_BEGIN name" and "IMG_END name".
shot() {
    screencapture -x "$OUT/$1.png" || { echo "screencapture failed for $1"; return; }
    if sips -s format jpeg -s formatOptions 60 -Z 1100 "$OUT/$1.png" --out "$OUT/$1.jpg" > /dev/null 2>&1; then
        echo "IMG_BEGIN $1"
        base64 -i "$OUT/$1.jpg" | tr -d '\n'
        echo
        echo "IMG_END $1"
    fi
}

dump_tab() {   # $1 = file prefix, $2 = element to click ("" = none)
    echo "::group::UI text: $1"
    with_timeout 40 osascript "$OUT/tab.applescript" "$2" | tee "$OUT/$1.txt"
    echo "::endgroup::"
    shot "$1"
}

echo "== launching $APP"
open -n "$APP"
sleep 15
pgrep -x "ΨLab" >/dev/null && echo "ΨLab process is running" || echo "ΨLab process NOT running"

shot "00-launch"
echo "::group::accessibility tree"
with_timeout 40 osascript "$OUT/tab.applescript" "__TREE__" | tee "$OUT/tree.txt"
echo "::endgroup::"
dump_tab "01-summary" ""
dump_tab "02-hartree-fock" "Hartree–Fock"
dump_tab "03-full-ci" "Full CI"
dump_tab "04-potential-curve" "Potential curve"

echo "== running the potential-curve scan"
echo "::group::accessibility tree (potential curve tab)"
with_timeout 40 osascript "$OUT/tab.applescript" "__TREE__" | tee "$OUT/tree-curve.txt"
echo "::endgroup::"
echo "::group::click Run scan"
with_timeout 90 osascript "$OUT/tab.applescript" "__SCAN__" 
echo "::endgroup::"
sleep 8
dump_tab "05-curve-after-scan" ""

echo "== menu commands"
echo "::group::menu items"
cat > "$OUT/menus.applescript" <<'APPLESCRIPT'
tell application "System Events"
    tell process "ΨLab"
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
with_timeout 30 osascript "$OUT/menus.applescript" | tee "$OUT/menus.txt"
echo "::endgroup::"

echo "== Edit > Copy Summary as Text"
pbcopy < /dev/null
cat > "$OUT/copy.applescript" <<'APPLESCRIPT'
tell application "System Events"
    tell process "ΨLab"
        click menu item "Copy Summary as Text" of menu 1 of menu bar item "Edit" of menu bar 1
    end tell
end tell
APPLESCRIPT
with_timeout 30 osascript "$OUT/copy.applescript" || echo "could not click the menu item"
sleep 1
echo "::group::clipboard after Copy Summary as Text"
pbpaste | tee "$OUT/clipboard-summary.txt"
echo "::endgroup::"

echo "== File > Export Results as Text (accept the save panel with Return)"
cat > "$OUT/export-text.applescript" <<'APPLESCRIPT'
tell application "System Events"
    tell process "ΨLab"
        set frontmost to true
        click menu item "Export Results as Text…" of menu 1 of menu bar item "File" of menu bar 1
        delay 3
        key code 36
        delay 3
    end tell
end tell
APPLESCRIPT
with_timeout 40 osascript "$OUT/export-text.applescript" || echo "could not drive the export"
EXPORTED=$(find "$HOME" -name 'PsiLab-HeH+-STO3G-R1.4632.txt' -mmin -10 2>/dev/null | head -n 1)
echo "::group::exported text file"
if [ -n "$EXPORTED" ]; then
    echo "file: $EXPORTED ($(wc -l < "$EXPORTED") lines)"
    cp "$EXPORTED" "$OUT/exported-results.txt"
    head -n 40 "$EXPORTED"
    echo "[...]"
    tail -n 60 "$EXPORTED"
else
    echo "exported file NOT found"
fi
echo "::endgroup::"

# --- dark mode pass: relaunch with the in-app appearance forced to dark ---------------
echo "== dark mode pass"
pkill -x "ΨLab" || true
sleep 3
open -n "$APP" --args -appearanceChoice dark
sleep 15
click_shot() {   # $1 = file prefix, $2 = tab to click ("" = none)
    with_timeout 60 osascript "$OUT/tab.applescript" "$2" > "$OUT/$1.txt" 2>&1
    shot "$1"
}
click_shot "dark-01-summary" "Summary"
click_shot "dark-02-hartree-fock" "Hartree–Fock"
click_shot "dark-03-full-ci" "Full CI"
click_shot "dark-04-potential-curve" "Potential curve"
with_timeout 90 osascript "$OUT/tab.applescript" "__SCAN__" | grep -E "CLICKED|NOT FOUND"
sleep 8
shot "dark-05-curve-after-scan"

# --- ΨLab 3.0: H2, LiH, H2O (launch arguments select molecule, basis set and method) ---------
# Selects a tab through View > <item> (unambiguous: the sidebar also has a "Full CI" segment),
# then dumps the visible text and takes a screenshot.
cat > "$OUT/view-menu.applescript" <<'APPLESCRIPT'
on run argv
    tell application "System Events"
        tell process "ΨLab"
            set frontmost to true
            click menu item (item 1 of argv) of menu 1 of menu bar item "View" of menu bar 1
        end tell
    end tell
    return "VIEW MENU " & (item 1 of argv)
end run
APPLESCRIPT
view_tab() {   # $1 = file prefix, $2 = View menu item
    with_timeout 20 osascript "$OUT/view-menu.applescript" "$2" || echo "could not select $2"
    sleep 2
    dump_tab "$1" ""
}
relaunch() {   # $@ = launch arguments
    pkill -x "ΨLab" || true
    sleep 3
    open -n "$APP" --args "$@"
    sleep 20
}
echo "== H2O, RHF/aug-cc-pVDZ"
relaunch -appearanceChoice light -psilabMolecule H2O -psilabBasis aug-cc-pVDZ -psilabMethod RHF
view_tab "10-h2o-rhf-summary" "Summary"
view_tab "11-h2o-rhf-hartree-fock" "Hartree–Fock"
view_tab "12-h2o-rhf-report" "Report"

echo "== LiH, Full CI/aug-cc-pVDZ"
relaunch -appearanceChoice light -psilabMolecule LiH -psilabBasis aug-cc-pVDZ -psilabMethod FCI
view_tab "13-lih-fci-summary" "Summary"
view_tab "14-lih-fci-full-ci" "Full CI"

echo "== H2, Full CI/aug-cc-pVDZ potential curve"
relaunch -appearanceChoice light -psilabMolecule H2 -psilabBasis aug-cc-pVDZ -psilabMethod FCI
view_tab "15-h2-curve" "Potential Curve"
with_timeout 90 osascript "$OUT/tab.applescript" "__SCAN__" | grep -E "CLICKED|NOT FOUND"
sleep 25
dump_tab "16-h2-curve-after-scan" ""

echo "== H2O, Full CI (must be refused) and a long calculation (warning)"
relaunch -appearanceChoice light -psilabMolecule H2O -psilabBasis aug-cc-pVDZ -psilabMethod FCI
view_tab "17-h2o-fci-refused" "Summary"
relaunch -appearanceChoice dark -psilabMolecule H2O -psilabBasis aug-cc-pCV5Z -psilabMethod RHF
view_tab "18-h2o-cv5z-warning-dark" "Summary"

# --- ΨLab 4.0: coupled cluster ------------------------------------------------------------
cat > "$OUT/method-menu.applescript" <<'APPLESCRIPT'
on run argv
    tell application "System Events"
        tell process "ΨLab"
            set frontmost to true
            set els to entire contents of window 1
            repeat with e in els
                try
                    if (role of e) is "AXPopUpButton" and (name of e) is "Method" then
                        click e
                        return "CLICKED Method pop-up"
                    end if
                end try
            end repeat
        end tell
    end tell
    return "NOT FOUND Method pop-up"
end run
APPLESCRIPT
echo "== H2O, CCSD(T)/aug-cc-pVDZ and the method menu"
relaunch -appearanceChoice light -psilabMolecule H2O -psilabBasis aug-cc-pVDZ -psilabMethod "CCSD(T)"
view_tab "19-h2o-ccsd-t-summary" "Summary"
view_tab "20-h2o-ccsd-t-coupled-cluster" "Coupled Cluster"
with_timeout 30 osascript "$OUT/method-menu.applescript" | grep -E "CLICKED|NOT FOUND"
sleep 2
shot "21-method-menu"
osascript -e 'tell application "System Events" to key code 53' || true

echo "== LiH, CCSDT/aug-cc-pVDZ (dark) and H2O CCSDTQ/aug-cc-pVTZ (refused)"
relaunch -appearanceChoice dark -psilabMolecule LiH -psilabBasis aug-cc-pVDZ -psilabMethod CCSDT
view_tab "22-lih-ccsdt-coupled-cluster-dark" "Coupled Cluster"
relaunch -appearanceChoice light -psilabMolecule H2O -psilabBasis aug-cc-pVTZ -psilabMethod CCSDTQ
view_tab "23-h2o-ccsdtq-refused" "Summary"

pkill -x "ΨLab" || true

# Crash reports of this run (none expected): print the crashed thread's backtrace.
echo "::group::crash reports"
find "$HOME/Library/Logs/DiagnosticReports" /Library/Logs/DiagnosticReports -name '*Lab*' -mmin -40 2>/dev/null | while read -r report; do
    echo "== $report"
    python3 - "$report" <<'PY' || head -c 6000 "$report"
import json, sys
text = open(sys.argv[1]).read()
header, _, body = text.partition("\n")
try:
    d = json.loads(body)
except Exception:
    print(text[:6000]); sys.exit(0)
print("exception:", d.get("exception"), "termination:", d.get("termination", {}).get("indicator"))
print("asi:", d.get("asi"))
images = d.get("usedImages", [])
for t in d.get("threads", []):
    if t.get("triggered"):
        for f in t.get("frames", [])[:40]:
            img = images[f.get("imageIndex", 0)].get("name", "?") if images else "?"
            print("  ", img, f.get("symbol", "?"), "+", f.get("symbolLocation", ""))
PY
done
echo "::endgroup::"
exit 0
