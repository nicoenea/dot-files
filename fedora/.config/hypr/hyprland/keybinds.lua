---------------------
---- KEYBINDINGS ----
---------------------
-- exec_cmd runs through `sh -c`, so $HOME expands inside command strings.
local mod      = "SUPER"
local terminal = "ghostty"
-- caelestia-shell has a built-in launcher — toggle it via the CLI.
-- Fallback: wofi --show drun
local launcher = "$HOME/.local/bin/caelestia shell drawers toggle launcher"

local exec = hl.dsp.exec_cmd

-- Core
hl.bind(mod .. " + Return",           exec(terminal))
hl.bind(mod .. " + Space",            exec(launcher))
hl.bind(mod .. " + C",                hl.dsp.window.close())
hl.bind(mod .. " + SHIFT + CTRL + Q", hl.dsp.exit())
hl.bind(mod .. " + F",                hl.dsp.window.fullscreen({ mode = "fullscreen" }))
hl.bind(mod .. " + SHIFT + F",        hl.dsp.window.fullscreen({ mode = "maximized" }))
hl.bind(mod .. " + V",                hl.dsp.window.float({ action = "toggle" }))
hl.bind(mod .. " + P",                hl.dsp.window.pseudo())
hl.bind(mod .. " + E",                hl.dsp.layout("togglesplit")) -- dwindle only

-- OpenRGB — toggle RGB lighting control window
hl.bind(mod .. " + O", hl.dsp.workspace.toggle_special("openrgb"))

-- Game mode — toggles caelestia game mode (suppresses notifs) and hides/shows the bar
hl.bind(mod .. " + G", exec("$HOME/.local/bin/caelestia-gamemode"))

-- Notification center — caelestia's sidebar drawer (popups expire after 5s; history lives here)
hl.bind(mod .. " + N", hl.dsp.global("caelestia:sidebar"))

-- Rapid-fire toggle: hold the button to rapid-fire while active (~20 clicks/sec)
hl.bind(mod .. " + R", exec("$HOME/.local/bin/rapidfire-toggle left"))
hl.bind(mod .. " + M", exec("$HOME/.local/bin/rapidfire-toggle middle"))

-- Screenshot (macOS-style)
hl.bind("CTRL + SHIFT + 3", exec("grim - | wl-copy"))
hl.bind("CTRL + SHIFT + 4", exec('grim -g "$(slurp)" - | wl-copy'))

-- Clipboard history
hl.bind(mod .. " + semicolon", exec("cliphist list | wofi --dmenu | cliphist decode | wl-copy"))

-- Focus (ALT + vim keys / arrows) and move windows (ALT + SHIFT + same keys)
local directions = {
    { "h", "left", "left" }, { "l", "right", "right" },
    { "k", "up",   "up" },   { "j", "down",  "down" },
}
for _, d in ipairs(directions) do
    local vim, arrow, dir = d[1], d[2], d[3]
    hl.bind("ALT + " .. vim,           hl.dsp.focus({ direction = dir }))
    hl.bind("ALT + " .. arrow,         hl.dsp.focus({ direction = dir }))
    hl.bind("ALT + SHIFT + " .. vim,   hl.dsp.window.move({ direction = dir }))
    hl.bind("ALT + SHIFT + " .. arrow, hl.dsp.window.move({ direction = dir }))
end

-- Resize (arrow keys, repeats while held)
hl.bind(mod .. " + CTRL + right", hl.dsp.window.resize({ x = 40,  y = 0,   relative = true }), { repeating = true })
hl.bind(mod .. " + CTRL + left",  hl.dsp.window.resize({ x = -40, y = 0,   relative = true }), { repeating = true })
hl.bind(mod .. " + CTRL + up",    hl.dsp.window.resize({ x = 0,   y = -40, relative = true }), { repeating = true })
hl.bind(mod .. " + CTRL + down",  hl.dsp.window.resize({ x = 0,   y = 40,  relative = true }), { repeating = true })

-- Workspaces: ASUS (main, left) on ALT+1-5, Gigabyte (secondary, right) on ALT+Q-T.
-- SHIFT+ALT+same key moves the active window there.
local ws_keys = { "1", "2", "3", "4", "5", "Q", "W", "E", "R", "T" }
for ws, key in ipairs(ws_keys) do
    hl.bind("ALT + " .. key,         hl.dsp.focus({ workspace = ws }))
    hl.bind("SHIFT + ALT + " .. key, hl.dsp.window.move({ workspace = ws }))
end

-- Scroll through workspaces
hl.bind(mod .. " + mouse_down", hl.dsp.focus({ workspace = "e-1" }))
hl.bind(mod .. " + mouse_up",   hl.dsp.focus({ workspace = "e+1" }))

-- Move/resize with mouse
hl.bind(mod .. " + mouse:272", hl.dsp.window.drag(),   { mouse = true })
hl.bind(mod .. " + mouse:273", hl.dsp.window.resize(), { mouse = true })

-- Move focus between monitors
hl.bind(mod .. " + comma",  hl.dsp.focus({ monitor = "l" }))
hl.bind(mod .. " + period", hl.dsp.focus({ monitor = "r" }))

-- Move workspace to other monitor
hl.bind(mod .. " + SHIFT + comma",  hl.dsp.workspace.move({ monitor = "l" }))
hl.bind(mod .. " + SHIFT + period", hl.dsp.workspace.move({ monitor = "r" }))

-- Media keys
hl.bind("XF86AudioRaiseVolume",  exec("wpctl set-volume @DEFAULT_AUDIO_SINK@ 5%+"))
hl.bind("XF86AudioLowerVolume",  exec("wpctl set-volume @DEFAULT_AUDIO_SINK@ 5%-"))
hl.bind("XF86AudioMute",         exec("wpctl set-mute @DEFAULT_AUDIO_SINK@ toggle"))
hl.bind("XF86AudioPlay",         exec("playerctl play-pause"))
hl.bind("XF86AudioNext",         exec("playerctl next"))
hl.bind("XF86AudioPrev",         exec("playerctl previous"))
hl.bind("XF86MonBrightnessUp",   exec("brightnessctl set 5%+"))
hl.bind("XF86MonBrightnessDown", exec("brightnessctl set 5%-"))

-- Lock screen
hl.bind(mod .. " + Escape", exec("caelestia shell lock lock"))
hl.bind(mod .. " + L",      exec("caelestia shell lock lock"))
