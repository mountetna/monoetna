export const RUN_ONCE = 0;
export const RUN_NEVER = -1;
export const RUN_INTERVAL = 1;

export const getRunState = (run) => (run != RUN_ONCE && run != RUN_NEVER ? RUN_INTERVAL : run);
export const getRunIntervalTime = (run) => (run != RUN_ONCE && run != RUN_NEVER ? run : 600);

export const formatTime = (time) => Intl.DateTimeFormat('en-US', {
   year: 'numeric', month: 'short', day: 'numeric',
   hour: 'numeric', minute: 'numeric'
}).format(new Date(time))

export const runTime = (ran_at, run) => {
  if (run == RUN_NEVER) return 'never';
  if (run == RUN_ONCE || !ran_at) return 'pending';

  const ranAt = new Date(ran_at);
  const runAt = new Date( ranAt.getTime() + run * 1000 );
  const now = new Date();

  if (runAt < now) return 'pending';

  return formatTime( runAt );
}

