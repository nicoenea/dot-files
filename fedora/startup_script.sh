# Visual Studio Code
sudo rpm --import https://packages.microsoft.com/keys/microsoft.asc
sudo sh -c 'echo -e "[code]\nname=Visual Studio Code\nbaseurl=https://packages.microsoft.com/yumrepos/vscode\nenabled=1\ngpgcheck=1\ngpgkey=https://packages.microsoft.com/keys/microsoft.asc" > /etc/yum.repos.d/vscode.repo'
dnf check-update
sudo dnf install code # or code-insiders

# pyenv
sudo yum -y install git gcc zlib-devel bzip2-devel readline-devel sqlite-devel openssl-devel
git clone https://github.com/pyenv/pyenv.git $HOME/.pyenv
# Zsh
echo 'Installing ZSH'
sudo dnf install zsh

echo 'Making ZSH as default shell...'
echo '\n Current shell setup: '
grep $USER /etc/passwd
chsh -s $(which zsh)
echo '\n new shell setup: '
grep $USER /etc/passwd

# Powerlevel10k

git clone --depth=1 https://github.com/romkatv/powerlevel10k.git ${ZSH_CUSTOM:-$HOME/.oh-my-zsh/custom}/themes/powerlevel10k
# install recommended font https://github.com/romkatv/powerlevel10k?tab=readme-ov-file#meslo-nerd-font-patched-for-powerlevel10k

# Oh-My-ZSH
echo 'Installing Oh-My-Zsh'
sh -c "$(curl -fsSL https://raw.github.com/ohmyzsh/ohmyzsh/master/tools/install.sh)"

#ohmyzsh plugins
git clone https://github.com/zsh-users/zsh-autosuggestions ${ZSH_CUSTOM:-~/.oh-my-zsh/custom}/plugins/zsh-autosuggestions
git clone https://github.com/zsh-users/zsh-syntax-highlighting.git ${ZSH_CUSTOM:-~/.oh-my-zsh/custom}/plugins/zsh-syntax-highlighting

