type Maybe<T> = undefined | T
type SMap<T> = {[k: string]: T}

type Nominal<S, T> = {_____nominal: T, s: S};
type KV<op> = Nominal<['kv', op], string>;
type Label = KV<'='>
type Constraint = KV<'=='>
type Duration = Nominal<'duration', string>;

export function duration(num: number, period: 's' = 's'): Duration {
  return `${Math.floor(num)}${period}` as any;
}

function kv<Op extends string>(k: string, v: string, op: Op): KV<Op> {
  if (k.includes('=') || v.includes('=')) {
    throw new Error(`Invalid key or value: ${k} -- ${v}`)
  }

  return `${k}${op}${v}` as any;
}

export function label(k: string, v: string): Label {
  return kv(k, v, '=')
}

export function constraint(k: string, v: string): Constraint {
  return kv(k, v, '==')
}

export const defaultSecretMeta = {
  external: undefined as Maybe<boolean>,
  name: undefined as Maybe<string>,
}

export type SecretMeta = typeof defaultSecretMeta;

export const defaultNetworkMeta = {
  external: undefined as Maybe<boolean>,
  name: undefined as Maybe<string>,
}

export type NetworkMeta = typeof defaultNetworkMeta;

export const defaultConfigMeta = {
  external: undefined as Maybe<boolean>,
}

export type ConfigMeta = typeof defaultConfigMeta;

export const defaultConfig = {
  source: '',
  target: '',
}

export type Config = typeof defaultConfig;

export const defaultNetwork = {
  aliases: undefined as Maybe<string[]>,
}

export type Network = typeof defaultNetwork;

export const defaultSecret = {
  source: '',
  target: undefined as Maybe<string>,
  uid: undefined as Maybe<string>,
  gid: undefined as Maybe<string>,
  mode: undefined as Maybe<string>,
}

export type Secret = typeof defaultSecret;

export const defaultDeploy = {
  update_config: {
    parallelism: 1,
    order: 'stop-first' as 'start-first' | 'stop-first',
  },
  labels: [] as Label[],
  placement: {
    constraints: [] as Constraint[],
  },
  configs: [] as Config[],
  networks: {} as SMap<Network>,
  secrets: {} as SMap<Secret>,
}

export const defaultHealthCheck = {
  interval: undefined as Maybe<string>,
  start_period: undefined as Maybe<string>,
  retries: 3,
  test: [] as string[],
  timeout: duration(3),
};

export type HealthCheck = typeof defaultHealthCheck;

export const defaultVolume = {
  type: 'volume' as 'volume' | 'bind',
  source: '',
  target: '',
  volume: {
    nocopy: false,
  }
}

export type Volume = typeof defaultVolume;

export const defaultService = {
  environment: {} as SMap<string>,
  image: "",
  volumes: [] as Volume[],
  user: undefined as Maybe<String>,
  command: [] as string[],
  healthcheck: undefined as Maybe<HealthCheck>
}

export type Service = typeof defaultService;

export const defaultStack = {
  version: '3.8',
  services: {} as SMap<Service>,
  secrets: {} as SMap<SecretMeta>,
  configs: {} as SMap<ConfigMeta>,
  networks: {} as SMap<NetworkMeta>,
}