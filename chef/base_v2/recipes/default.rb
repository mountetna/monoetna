#
# Cookbook:: base_v2
# Recipe:: default
#
# Copyright:: 2020, The Authors, All Rights Reserved.
include_recipe 'base_v2::packages'
include_recipe 'base_v2::pgsql'
include_recipe 'base_v2::users'
include_recipe 'base_v2::docker'
