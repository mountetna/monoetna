load 'bin/metis_client'

describe MetisShell do
  describe MetisShell::Ls do
    before(:each) do
      MetisConfig.instance.config("tmp-config")
    end
    it 'lists a directory' do
      expect{MetisShell.new("metis://metis.test/labors", "ls").run}.to output(/No prokect is selected/).to_stdout
    end

    after(:each) do
      ::File.unlink("tmp-config")
    end
  end
end
