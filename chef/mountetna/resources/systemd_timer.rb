property :unit, String, required: true
property :time_seconds, Integer, required: true

action :create do
  systemd_unit "#{new_resource.name}.timer" do
    content({
        Unit: { Description: "Starts the #{new_resource.unit} job every #{new_resource.time_seconds} seconds" },
        Timer: { OnActiveSec: new_resource.time_seconds, Unit: new_resource.unit, "OnUnitActiveSec": new_resource.time_seconds, },
        Install: { WantedBy: "multi-user.target" },
    })

    action [:create, :enable, :restart]
  end
end
