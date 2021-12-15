describe Metis::ThumbnailCache do
  before(:each) do
    stubs.clear("data_blocks")
    bucket = default_bucket("athena")

    blueprints_folder = create_folder("athena", "blueprints")
    stubs.create_folder("athena", "files", "blueprints")

    @helmet_file = create_file("athena", "helmet.txt", HELMET, folder: blueprints_folder)
    stubs.create_file("athena", "files", "blueprints/helmet.txt", HELMET)

    red_square_png = "\x89PNG\r\n\u001A\n\u0000\u0000\u0000\rIHDR\u0000\u0000\u0000\u0004\u0000\u0000\u0000\u0004\b\u0002\u0000\u0000\u0000&\x93\t)\u0000\u0000\u0001\x84iCCPICC profile\u0000\u0000(\x91}\x91=H\xC3@\u001C\xC5_ӊE*\x82\xED \xE2\u0010\xA4:Y\u0010\u0015q\xD4*\u0014\xA1B\xA9\u0015Zu0\xB9\xF4\v\x9A4$).\x8E\x82k\xC1\xC1\x8FŪ\x83\x8B\xB3\xAE\u000E\xAE\x82 \xF8\u0001\xE2\xE4\xE8\xA4\xE8\"%\xFE/)\xB4\x88\xF1\xE0\xB8\u001F\xEF\xEE=\xEE\xDE\u0001B\xA3\xC2T30\u000E\xA8\x9Ae\xA4\u0013q1\x9B[\u0015\xBB_\u0011D\u0000a\xF4cXb\xA6>\x97J%\xE19\xBE\xEE\xE1\xE3\xEB]\x8Cgy\x9F\xFBs\xF4*y\x93\u0001>\x91x\x96\xE9\x86E\xBCA<\xBDi\xE9\x9C\xF7\x89#\xAC$)\xC4\xE7\xC4c\u0006]\x90\xF8\x91\xEB\xB2\xCBo\x9C\x8B\u000E\v<3bd\xD2\xF3\xC4\u0011b\xB1\xD8\xC1r\a\xB3\x92\xA1\u0012O\u0011G\u0015U\xA3|!\xEB\xB2\xC2y\x8B\xB3Z\xA9\xB1\xD6=\xF9\vCyme\x99\xEB4\x87\x90\xC0\"\x96\x90\x82\b\u00195\x94Q\x81\x85\u0018\xAD\u001A)&Ҵ\u001F\xF7\xF0\u000F:\xFE\u0014\xB9dr\x95\xC1ȱ\x80*TH\x8E\u001F\xFC\u000F~wk\u0016&'ܤP\u001C\xE8z\xB1\xED\x8F\u0011\xA0{\u0017h\xD6m\xFB\xFBض\x9B'\x80\xFF\u0019\xB8\xD2\xDA\xFEj\u0003\x98\xF9$\xBD\xDE֢G@\xDF6pq\xDD\xD6\xE4=\xE0r\a\u0018x\xD2%Cr$?M\xA1P\u0000\xDE\xCF\xE8\x9Br@\xF8\u0016\xE8Ys{k\xED\xE3\xF4\u0001\xC8PW\xC9\e\xE0\xE0\u0010\u0018-R\xF6\xBAǻ\x83\x9D\xBD\xFD{\xA6\xD5\xDF\u000F9nr\x90݁ږ\u0000\u0000\u0000\tpHYs\u0000\u0000.#\u0000\u0000.#\u0001x\xA5?v\u0000\u0000\u0000\atIME\a\xE4\b\u0004\u0013\u0014&\u001C\xD0`#\u0000\u0000\u0000\u0014IDAT\b\xD7c|\xAE\xAA\xCA\u0000\u0003L\fH\u00007\a\u0000=X\u00019\u0014\u001Dh\u001C\u0000\u0000\u0000\u0000IEND\xAEB`\x82"
    stubs.create_file("athena", "files", "red_square.png", red_square_png)
    @target_file = create_file("athena", "red_square.png", red_square_png, bucket: bucket)

    @cache = Metis::ThumbnailCache.new
  end

  after(:each) do
    stubs.clear("data_blocks")

    expect(stubs.contents(:athena)).to be_empty
  end

  it "correctly identifies image file" do
    expect(@cache.image?(@helmet_file)).to eq(false)
    expect(@cache.image?(@target_file)).to eq(true)
  end

  it "can generate a thumbnail" do
    expect(@cache.thumbnail_in_cache?(@target_file)).to eq(false)

    @cache.generate_thumbnail(@target_file)

    expect(@cache.thumbnail_in_cache?(@target_file)).to eq(true)
  end

  it "can fetch a cached thumbnail" do
    expect {
      @cache.thumbnail(@target_file)
    }.to raise_error(Metis::ThumbnailNotExistError)

    @cache.generate_thumbnail(@target_file)

    expect {
      @cache.thumbnail(@target_file)
    }.not_to raise_error(Metis::ThumbnailError)
  end
end
