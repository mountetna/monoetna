Sequel.migration do
  up do
    fetch('SELECT MIN(id) as id, folder_id, folder_name, bucket_id, count(*) as cnt FROM folders GROUP BY folder_id, folder_name, bucket_id HAVING count(id) > 1').each do |f|
      folder_ids = from(:folders).where(folder_id: f[:folder_id], folder_name: f[:folder_name], bucket_id: f[:bucket_id]).map { |r| r[:id] }

      from(:folders).where(folder_id: folder_ids).update(folder_id: f[:id])
      from(:files).where(folder_id: folder_ids).update(folder_id: f[:id])
      from(:folders).where(id: folder_ids.select { |id| id != f[:id] }).delete
    end
  end
end
