/** @type {import('next').NextConfig} */
const nextConfig = {
    output: "standalone",
    // TODO: remove for prod
    images: {
        remotePatterns: [
            {
                protocol: 'https',
                hostname: '**'
            },
        ],
    },
};

export default nextConfig;
