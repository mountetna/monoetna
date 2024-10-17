/** @type {import('next').NextConfig} */
const nextConfig = {
    output: "standalone",
    images: {
        remotePatterns: [
            {
                protocol: 'https',
                hostname: '**'
            },
        ],
    },
    experimental: {
        serverActions: {
            allowedOrigins: [
                // this is for local development
                // TODO: see if this is ok for production
                // (i.e. for vesta.ucsf.edu, vesta.mountetna.dev, etc.)
                'vesta.development.local',
                'vesta_app_fe',
            ],
        }
    }
};

export default nextConfig;
