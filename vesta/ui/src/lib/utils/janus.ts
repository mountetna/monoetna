export interface Permission {
    role: string
    privileged: boolean
    project_name: string
}

const ROLES: Record<string, string> = { a: 'administrator', e: 'editor', v: 'viewer', g: 'guest' };

export const isSuperuser = ({ permissions }: { permissions: Record<string, Permission> }) => (permissions.administration && permissions.administration.role == ROLES.a)
export const isSuperEditor = ({ permissions }: { permissions: Record<string, Permission> }) => (permissions.administration && [ROLES.a, ROLES.e].includes(permissions.administration.role))
export const isSuperViewer = ({ permissions }: { permissions: Record<string, Permission> }) => (permissions.administration && [ROLES.a, ROLES.e, ROLES.v].includes(permissions.administration.role))
export const isEditor = ({ permissions }: { permissions: Record<string, Permission> }, project_name: string) => (permissions[project_name] && [ROLES.a, ROLES.e].includes(permissions[project_name].role)) || isSuperEditor({ permissions })
export const isAdmin = ({ permissions }: { permissions: Record<string, Permission> }, project_name: string) => (permissions[project_name] && permissions[project_name].role == ROLES.a) || isSuperuser({ permissions })
export const isPrivileged = ({ permissions }: { permissions: Record<string, Permission> }, project_name: string) => (permissions[project_name].privileged);
export const isGuest = ({ permissions }: { permissions: Record<string, Permission> }, project_name: string) => (permissions[project_name] && permissions[project_name].role == ROLES.g)

function parsePermissions(perms: string) {
    // Permissions are encoded as 'a:project1,project2;v:project3'
    return perms.split(/;/).map(perm => {
        if (!perm) return [];
        let [encoded_role, projects] = perm.split(/:/);
        let role = ROLES[encoded_role.toLowerCase()];
        let privileged = encoded_role.toUpperCase() == encoded_role;
        return projects.split(/,/).map(
            project_name => ({ role, privileged, project_name })
        )
    }).reduce((perms, projects) => {
        projects.forEach(proj => perms[proj.project_name] = proj);
        return perms
    }, ({} as Record<string, Permission>));
}

const parseFlags = (flags: string) => flags ? flags.split(';') : []

export function parseToken(token: string) {
    let [header, params, signature] = token.split(/\./);
    let {
        email,
        name,
        perm,
        flags,
        joined_at,
    }: {
        email: string,
        name: string,
        perm: string,
        flags: string,
        joined_at: string,
    } = JSON.parse(atob(params));

    return {
        email,
        name,
        permissions: parsePermissions(perm),
        flags: parseFlags(flags),
        joined_at: new Date(joined_at),
    };
}
