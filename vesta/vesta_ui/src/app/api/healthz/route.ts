export async function GET(request: Request) {
    return new Response('Application is healthy', {
        status: 200,
    })
}