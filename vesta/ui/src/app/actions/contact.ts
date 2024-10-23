'use server'

import { z } from 'zod'

import { VestaApiClient } from "@/lib/clients/vesta-api/client"
import { SendContactStatus } from "@/lib/clients/vesta-api/models"
import rateLimit from '@/lib/utils/rate-limit'


const emailFormSchema = z.object({
    email: z.string({
        invalid_type_error: 'Invalid Email',
    }),
})

const rateLimitMax = parseInt(process.env.CONTACT_EMAIL_RATE_LIMIT)
const rateLimitIntervalSeconds = parseInt(process.env.CONTACT_EMAIL_RATE_LIMIT_INTERVAL_SECONDS)

const rateLimiter = rateLimit({
    interval: rateLimitIntervalSeconds * 1000, // convert from milliseconds
    uniqueTokenPerInterval: 1,
})

const rateLimiterToken = 'RATE_LIMIT'


export async function sendContributeEmail(formData: FormData): Promise<SendContactStatus> {
    const validatedFields = emailFormSchema.safeParse({
        email: formData.get('email'),
    })

    if (!validatedFields.success) {
        return {
            status: 'error',
            message: String(validatedFields.error.flatten().fieldErrors.email),
        }
    }

    try {
        await rateLimiter.check(rateLimitMax, rateLimiterToken)
    } catch {
        return {
            status: 'rateLimit',
        }
    }

    const apiClient = new VestaApiClient()
    return await apiClient.sendContributeEmail(validatedFields.data.email)
}