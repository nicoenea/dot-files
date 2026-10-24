-- Hyprland >= 0.55 Lua config. Modules live in ./hyprland/*.lua (resolved
-- relative to this file). API reference: /usr/share/hypr/stubs/hl.meta.lua
-- and https://wiki.hypr.land

local ASUS     = "desc:ASUSTek COMPUTER INC VG27AC1A"
local GIGABYTE = "desc:GIGA-BYTE TECHNOLOGY CO. LTD. M27Q"

------------------
---- MONITORS ----
------------------
-- Matched by EDID description, not connector name. The DRM connector for the
-- Gigabyte shifted DP-3 -> DP-1 during the Aug 2026 RTX 3080 -> RX 9070 XT swap,
-- which silently dropped it to 60 Hz because no monitor= line matched. desc:
-- matching is tied to the panel's EDID, so it survives port and GPU changes.
hl.monitor({ output = ASUS,     mode = "2560x1440@170", position = "0x0",    scale = 1 }) -- main, left
hl.monitor({ output = GIGABYTE, mode = "2560x1440@170", position = "2560x0", scale = 1 }) -- secondary, right

-- Pin workspaces to monitors: 1-5 on the ASUS, 6-10 on the Gigabyte
for ws = 1, 10 do
    hl.workspace_rule({
        workspace = tostring(ws),
        monitor   = ws <= 5 and ASUS or GIGABYTE,
        default   = (ws == 1 or ws == 6) or nil,
    })
end

require("hyprland.env")
require("hyprland.execs")
require("hyprland.look")
require("hyprland.input")
require("hyprland.keybinds")
require("hyprland.rules")
