module Etna::Spec
  module Auth
    AUTH_USERS = {
      superuser: {
        email: 'zeus@olympus.org', first: 'Zeus', perm: 'A:administration'
      },
      admin: {
        email: 'hera@olympus.org', first: 'Hera', perm: 'a:labors'
      },
      editor: {
        email: 'eurystheus@twelve-labors.org', first: 'Eurystheus', perm: 'E:labors'
      },
      restricted_editor: {
        email: 'copreus@twelve-labors.org', first: 'Copreus', perm: 'e:labors'
      },
      viewer: {
        email: 'hercules@twelve-labors.org', first: 'Hercules', perm: 'v:labors'
      },
      non_user: {
        email: 'nessus@centaurs.org', first: 'Nessus', perm: ''
      }
    }

    def auth_header(user_type)
      header(*Etna::TestAuth.token_header(AUTH_USERS[user_type]))
    end
  end
end
