-----------------------
---- LOOK AND FEEL ----
-----------------------
hl.config({
    general = {
        gaps_in          = 4,
        gaps_out         = 8,
        border_size      = 2,
        col              = {
            active_border   = { colors = { "rgba(cba6f7ff)", "rgba(89b4faff)" }, angle = 45 },
            inactive_border = "rgba(313244ff)",
        },
        resize_on_border = true,
        layout           = "dwindle",
    },

    decoration = {
        rounding         = 8,
        active_opacity   = 1.0,
        inactive_opacity = 0.95,

        shadow = {
            enabled      = true,
            range        = 12,
            render_power = 3,
            color        = "rgba(1a1a2ecc)",
        },

        blur = {
            enabled           = true,
            size              = 6,
            passes            = 3,
            new_optimizations = true,
            xray              = false,
        },
    },

    animations = {
        enabled = true,
    },

    -- dwindle:pseudotile was removed in 0.55; the pseudo dispatcher (Super+P) still works.
    dwindle = {
        preserve_split = true,
    },

    master = {
        new_status = "slave",
    },

    misc = {
        force_default_wallpaper  = 0,
        disable_hyprland_logo    = true,
        disable_splash_rendering = true,
        mouse_move_enables_dpms  = true,
        key_press_enables_dpms   = true,

        -- Adaptive sync (FreeSync). Both panels report a 48-170 Hz EDID range, and
        -- AMD's VRR on Wayland is reliable in a way NVIDIA's was not. 2 = fullscreen
        -- apps only: full-time VRR (1) makes some panels flicker on the desktop,
        -- because the refresh rate swings wildly as the cursor moves.
        vrr = 2,
    },
})

hl.curve("easeOut", { type = "bezier", points = { { 0.16, 1 },    { 0.3, 1 } } })
hl.curve("easeIn",  { type = "bezier", points = { { 0.7, 0 },     { 0.84, 0 } } })
hl.curve("spring",  { type = "bezier", points = { { 0.25, 0.46 }, { 0.45, 0.94 } } })

hl.animation({ leaf = "windows",    enabled = true, speed = 4, bezier = "easeOut", style = "slide" })
hl.animation({ leaf = "windowsOut", enabled = true, speed = 3, bezier = "easeIn",  style = "slide" })
hl.animation({ leaf = "border",     enabled = true, speed = 6, bezier = "default" })
hl.animation({ leaf = "fade",       enabled = true, speed = 4, bezier = "easeOut" })
hl.animation({ leaf = "workspaces", enabled = true, speed = 4, bezier = "spring",  style = "slidevert" })
