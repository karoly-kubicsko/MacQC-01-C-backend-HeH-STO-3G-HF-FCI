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
    set logText to ""
    tell application "System Events"
        tell process "MacQC"
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
                -- that follows the "Step (bohr)" pop-up in the tree
                set els to entire contents of window 1
                set seenStep to false
                repeat with e in els
                    try
                        set r to role of e
                        set n to ""
                        try
                            set n to name of e
                        end try
                        if r is "AXPopUpButton" and n is "Step (bohr)" then
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

shot() { screencapture -x "$OUT/$1.png" || echo "screencapture failed for $1"; }

dump_tab() {   # $1 = file prefix, $2 = element to click ("" = none)
    echo "::group::UI text: $1"
    with_timeout 40 osascript "$OUT/tab.applescript" "$2" | tee "$OUT/$1.txt"
    echo "::endgroup::"
    shot "$1"
}

echo "== launching $APP"
open -n "$APP"
sleep 15
pgrep -x MacQC >/dev/null && echo "MacQC process is running" || echo "MacQC process NOT running"

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
with_timeout 30 osascript "$OUT/menus.applescript" | tee "$OUT/menus.txt"
echo "::endgroup::"

echo "== Edit > Copy Summary as Text"
pbcopy < /dev/null
cat > "$OUT/copy.applescript" <<'APPLESCRIPT'
tell application "System Events"
    tell process "MacQC"
        click menu item "Copy Summary as Text" of menu 1 of menu bar item "Edit" of menu bar 1
    end tell
end tell
APPLESCRIPT
with_timeout 30 osascript "$OUT/copy.applescript" || echo "could not click the menu item"
sleep 1
echo "::group::clipboard after Copy Summary as Text"
pbpaste | tee "$OUT/clipboard-summary.txt"
echo "::endgroup::"

# --- dark mode pass: relaunch with the in-app appearance forced to dark ---------------
echo "== dark mode pass"
pkill -x MacQC || true
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

# Print downscaled screenshots into the log as base64 JPEG so they can be read through the
# GitHub API without downloading the artifact. Decode: lines between IMG_BEGIN/IMG_END.
echo "== screenshots (base64 JPEG)"
for f in "$OUT"/*.png; do
    name=$(basename "$f" .png)
    sips -s format jpeg -s formatOptions 60 -Z 1100 "$f" --out "$OUT/$name.jpg" > /dev/null 2>&1 || continue
    echo "IMG_BEGIN $name"
    base64 -i "$OUT/$name.jpg" | tr -d '\n'
    echo
    echo "IMG_END $name"
done

pkill -x MacQC || true
exit 0
