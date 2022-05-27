class CcAgreement < Sequel::Model

  def to_hash
    {
      user_email: user_email,
      cc_text: cc_text,
      agreed: agreed,
      project_name: project_name
    }
  end
end
