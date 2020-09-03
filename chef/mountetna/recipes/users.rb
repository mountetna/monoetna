def run
  node['users'].each do |user_config|
    user_name = user_config[:name]

    user user_name do
      home "/home/#{user_name}"
      password user_config[:password]
      shell '/bin/bash'
      not_if "getent passwd #{user_name}"
    end

    directory "/home/#{user_name}" do
      action :create
      owner user_name
      group user_name
      mode '0700'
    end

    directory "/home/#{user_name}/.ssh" do
      action :create
      owner user_name
      group user_name
      mode '0700'
    end

    if user_config[:ssh_public_key]
      file "/home/#{user_name}/.ssh/authorized_keys" do
        content user_config[:ssh_public_key]
        mode '0600'
        owner user_name
        group user_name
        not_if { ::File.exists? "/home/#{user_name}/.ssh/authorized_keys" }
      end
    end

    file "/home/#{user_name}/README" do
      content user_config[:password]
      mode '0600'
      owner user_name
      group user_name
      not_if { ::File.exists? "/home/#{user_name}/README" }
    end

    group 'wheel' do
      action :modify
      append true
      members [user_name]
    end

    group 'docker'
    group 'docker' do
      action :modify
      append true
      members [user_name]
    end

    setup_dotfiles(user_name)
    setup_vim(user_name)
  end
end

def setup_dotfiles(usr)
  template("/home/#{usr}/.bashrc") do
    source('userbashrc.erb')
    owner(usr)
    group(usr)
    mode('0644')
  end

  template("/home/#{usr}/.bash_profile") do
    source('bashprofile.erb')
    owner(usr)
    group(usr)
    mode('0644')
  end

  if usr == 'root'
    template("/home/#{usr}/.tmux.conf") do
      source('tmux.conf.erb')
      owner(usr)
      group(usr)
      mode('0644')
    end
  end
end

def setup_vim(usr)
  template("/home/#{usr}/.vimrc") do
    source('vimrc.erb')
    owner(usr)
    group(usr)
    mode('0644')
  end

  setup_vim_folders(usr)

  git("/home/#{usr}/.vim/bundle/Vundle.vim") do
    repository('https://github.com/VundleVim/Vundle.vim.git')
    reference('master')
    action(:checkout)
  end
end

def setup_vim_folders(usr)
  directory("/home/#{usr}/.vim") do
    owner(usr)
    group(usr)
    mode('0755')
    action(:create)
  end

  directory("/home/#{usr}/.vim/tmp") do
    owner(usr)
    group(usr)
    mode('0755')
    action(:create)
  end

  directory("/home/#{usr}/.vim/bundle") do
    owner(usr)
    group(usr)
    mode('0755')
    action(:create)
  end
end

run
