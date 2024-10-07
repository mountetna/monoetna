require_relative 'controller'

class ContactController < Vesta::Controller
  def send_contact_email
    config = Vesta.instance.config(:vesta)[:contact_email]
    requester_email = @params[:requester_email]

    # send to DL team
    begin
      puts "Sending internal email"

      internal_config = config[:internal]

      send_email(
        'UCSF Data Library',
        internal_config[:address],
        internal_config[:subject],
        internal_config[:template] % {email: requester_email},
      )

      puts "Successfully sent internal email to #{internal_config[:address]}"
    rescue error => e
      puts "Error sending internal email to #{internal_config[:address]}: #{e}"
      raise
    end

    # send to DL team
    begin
      puts "Sending external email"

      external_config = config[:external]

      send_email(
        requester_email,
        requester_email,
        external_config[:subject],
        external_config[:template],
      )

      puts "Successfully sent external email to #{external_config[:address]}"
    rescue error => e
      puts "Error sending external email to #{external_config[:address]}: #{e}"
      raise
    end

    success_json({ success: 200 })
  end
end
