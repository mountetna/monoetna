module Etna::Spec
  module Auth
    AUTH_USERS = {
      superuser: {
        email: 'zeus@olympus.org', name: 'Zeus', perm: 'A:administration'
      },
      admin: {
        email: 'hera@olympus.org', name: 'Hera', perm: 'a:labors'
      },
      editor: {
        email: 'eurystheus@twelve-labors.org', name: 'Eurystheus', perm: 'E:labors'
      },
      restricted_editor: {
        email: 'copreus@twelve-labors.org', name: 'Copreus', perm: 'e:labors'
      },
      viewer: {
        email: 'hercules@twelve-labors.org', name: 'Hercules', perm: 'v:labors'
      },
      non_user: {
        email: 'nessus@centaurs.org', name: 'Nessus', perm: ''
      }
    }

    def auth_header(user_type)
      header(*Etna::TestAuth.token_header(AUTH_USERS[user_type]))
    end
  end
end
