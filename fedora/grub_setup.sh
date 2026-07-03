#!/usr/bin/env bash
set -e

# Apply preferred GRUB settings on Fedora.
# Removes any custom theme and configures hidden timeout for fast boot.
# Hold Esc/Shift at boot to show the menu (e.g. to pick Windows).
# Must be run as root: sudo ./grub_setup.sh

if [ "$EUID" -ne 0 ]; then
    echo "ERROR: Run this script with sudo."
    exit 1
fi

echo "Removing any installed custom GRUB themes..."
rm -rf /boot/grub2/themes/Elegant-*
rm -rf /boot/grub2/themes/catppuccin-*

echo "Restoring /etc/default/grub..."
if [ -f /etc/default/grub.bak ]; then
    cp /etc/default/grub.bak /etc/default/grub
    echo "  Restored from backup."
else
    sed -i 's|^GRUB_THEME=.*||' /etc/default/grub
    sed -i 's|^GRUB_BACKGROUND=.*||' /etc/default/grub
    sed -i 's|^GRUB_GFXMODE=.*||' /etc/default/grub
    sed -i 's|^GRUB_FONT=.*||' /etc/default/grub
    sed -i 's|^#GRUB_TERMINAL_OUTPUT=.*|GRUB_TERMINAL_OUTPUT="console"|' /etc/default/grub
fi

echo "Applying preferred settings..."

# Show menu briefly then auto-boot — reliable on UEFI, easy to catch for Windows
if grep -q "GRUB_TIMEOUT_STYLE=" /etc/default/grub; then
    sed -i 's|.*GRUB_TIMEOUT_STYLE=.*|GRUB_TIMEOUT_STYLE=menu|' /etc/default/grub
else
    echo 'GRUB_TIMEOUT_STYLE=menu' >> /etc/default/grub
fi

sed -i 's|.*GRUB_TIMEOUT=.*|GRUB_TIMEOUT=5|' /etc/default/grub

echo "Regenerating GRUB config..."
grub2-mkconfig -o /boot/grub2/grub.cfg

echo "Setting DNF installonly_limit=2 (keep current kernel + 1 fallback)..."
if grep -q "^installonly_limit=" /etc/dnf/dnf.conf; then
    sed -i 's/^installonly_limit=.*/installonly_limit=2/' /etc/dnf/dnf.conf
else
    echo 'installonly_limit=2' >> /etc/dnf/dnf.conf
fi

echo "Removing old kernels..."
old_kernels=$(dnf repoquery --installonly --latest-limit=-2 -q 2>/dev/null)
if [ -n "$old_kernels" ]; then
    echo "$old_kernels" | xargs dnf remove -y
else
    echo "  Nothing to remove."
fi

echo "Done. Reboot to apply."
