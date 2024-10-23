import { NextRequest, NextResponse } from 'next/server';

export const config = {
    matcher: [
        /*
 * Match all request paths except for the ones starting with:
 * - api (API routes)
 * - _next/static (static files)
 * - _next/image (image optimization files)
 * - favicon.ico (favicon file)
 */
        '/((?!api/healthz))',
    ],
}

export function middleware(req: NextRequest) {
    const basicAuthHeader = req.headers.get('authorization')
    const basicAuthUser = process.env.BASIC_AUTH_USER
    const basicAuthPassword = process.env.BASIC_AUTH_PASSWORD

    if (basicAuthUser && basicAuthPassword) {
        if (basicAuthHeader) {
            const authVal = basicAuthHeader.split(' ')[1]
            const [user, pwd] = atob(authVal).split(':')

            if (user === basicAuthUser && pwd === basicAuthPassword) {
                return NextResponse.next()
            }
        }

        const url = req.nextUrl
        url.pathname = '/api/basicauth'

        return NextResponse.rewrite(url)
    }
}