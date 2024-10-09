namespace NodeJS {
    interface ProcessEnv {
        BASIC_AUTH_USER?: string
        BASIC_AUTH_PASSWORD?: string
        API_URL: string
        JANUS_URL: string
        JANUS_TOKEN_COOKIE_NAME: string
        TIMUR_URL: string
        CONTACT_EMAIL_RATE_LIMIT: number
        CONTACT_EMAIL_RATE_LIMIT_INTERVAL_SECONDS: number
    }
}