export async function GET(request: Request) {
    return new Response('Authentication required', {
        status: 401,
        headers: {
            "WWW-Authenticate": "Basic realm='private_pages'",
        },
    })
}