export const projectTypeFor = project => project.resource ? project.requires_agreement ? 'community' : 'resource' : 'active';
