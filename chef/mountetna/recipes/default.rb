#
# Cookbook:: mountetna
# Recipe:: default
#
# Copyright:: 2020, The Authors, All Rights Reserved.
include_recipe 'mountetna::packages'
include_recipe 'mountetna::pgsql'
include_recipe 'mountetna::users'
include_recipe 'mountetna::docker'
include_recipe 'mountetna::edge_apache'

