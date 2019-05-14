require_relative '../../models/archimedes'

class ArchimedesController <  Etna::Controller
  def consignment
    manifests = @params[:queries]

    if @params[:manifest_ids]
      manifest_ids = @params[:manifest_ids].uniq

      manifests = Manifest.where(id: manifest_ids).all

      missing_manifest_ids = manifest_ids - manifests.map(&:id)

      unless missing_manifest_ids.empty?
        raise Etna::BadRequest, "Could not find manifests #{missing_manifest_ids.join(', ')}"
      end

      manifests = manifests.map(&:script)
    end

    raise 'No manifests requested!' unless manifests && !manifests.empty?

    consignments = manifests.map do |script|
      [
        # the identifier is the md5 hash of the script
        Digest::MD5.hexdigest(script),

        # the consignment is the values computed by Archimedes
        Archimedes::Manifest.new(
          token,
          @params[:project_name],
          script
        ).payload
      ]
    end.to_h

    success(consignments.to_json, 'application/json')
  rescue Archimedes::LanguageError => e
    raise Etna::BadRequest, e.message.to_s
  end

  private

  def with_record_name(script)
    @params[:record_name] ? "@record_name = '#{@params[:record_name].gsub(/'/, "''")}'\n#{script}" : script
  end

  def token
    @token ||= @request.cookies[Archimedes.instance.config(:token_name)]
  end
end

