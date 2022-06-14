class CcAgreement < Sequel::Model
  plugin :timestamps, update_on_create: true

  def to_hash
    {
      user_email: user_email,
      cc_text: cc_text,
      agreed: agreed,
      project_name: project_name
    }
  end
end
